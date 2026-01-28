// Double checking is on going with thesis wrting and reading 
// (this is becuase there were some errors occur 21/1/2026 in jacobian matrices)
// So far this is the most speed simulation I gained. ( Analatical + Augmented Lagrange + IKF)
// Need to check IKF again
// no torque is added
// Compile:
// g++ -std=c++11 -O3 Fourbar_Complete_Fixed.cpp -o fourbar_server
//
// ./fourbar_server
// ============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <algorithm>

const double PI = 3.14159265358979323846;

// ============ STRUCTURES ============
struct SystemParams {
    double g, L1, L2, L3, L4;
    double m1, m2, m3;
    double I1, I2, I3;
    double T0, win;
};

struct State {
    double q[9];    // Position: [Rx1, Ry1, θ1, Rx2, Ry2, θ2, Rx3, Ry3, θ3]
    double dq[9];   // Velocity
};

// ============ INDIRECT KALMAN FILTER STATE ============
struct IndirectEKFState {
    double delta_x[9];      // Error state: [delta_q, delta_dq] (9-dimensional)
    double P[9][9];         // Error covariance
    double Q[9][9];         // Process noise covariance
    double R[9][9];         // Measurement noise covariance
};

struct LogData {
    std::vector<double> time;
    std::vector<double> theta_matlab, omega_matlab;
    std::vector<double> theta_al, omega_al;
    std::vector<double> theta_corrected, omega_corrected;
    std::vector<double> cycle_time_us;
};

LogData logData;

// ============ BYTE SWAPPING FOR NETWORK ============
double swapDouble(double value) {
    uint64_t temp;
    std::memcpy(&temp, &value, sizeof(double));
    temp = ((temp & 0xFF00000000000000ULL) >> 56) |
           ((temp & 0x00FF000000000000ULL) >> 40) |
           ((temp & 0x0000FF0000000000ULL) >> 24) |
           ((temp & 0x000000FF00000000ULL) >> 8) |
           ((temp & 0x00000000FF000000ULL) << 8) |
           ((temp & 0x0000000000FF0000ULL) << 24) |
           ((temp & 0x000000000000FF00ULL) << 40) |
           ((temp & 0x00000000000000FFULL) << 56);
    double result;
    std::memcpy(&result, &temp, sizeof(double));
    return result;
}

// ============ MATRIX UTILITIES ============
void invert9x9(double A[9][9], double Ainv[9][9]) {
    double temp[9][9];
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            temp[i][j] = A[i][j];
            Ainv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    for (int col = 0; col < 9; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 9; ++row)
            if (std::fabs(temp[row][col]) > std::fabs(temp[pivot][col]))
                pivot = row;
        
        if (std::fabs(temp[col][col]) < 1e-14) {
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j)
                    Ainv[i][j] = (i == j) ? 1e14 : 0.0;
            return;
        }
        
        if (pivot != col) {
            for (int j = 0; j < 9; ++j) {
                std::swap(temp[col][j], temp[pivot][j]);
                std::swap(Ainv[col][j], Ainv[pivot][j]);
            }
        }
        
        double diag = temp[col][col];
        for (int j = 0; j < 9; ++j) {
            temp[col][j] /= diag;
            Ainv[col][j] /= diag;
        }
        
        for (int row = 0; row < 9; ++row) {
            if (row != col) {
                double factor = temp[row][col];
                for (int j = 0; j < 9; ++j) {
                    temp[row][j] -= factor * temp[col][j];
                    Ainv[row][j] -= factor * Ainv[col][j];
                }
            }
        }
    }
}

void mult9x9(double A[9][9], double B[9][9], double C[9][9]) {
    double temp[9][9] = {0};
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            for (int k = 0; k < 9; ++k)
                temp[i][j] += A[i][k] * B[k][j];
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            C[i][j] = temp[i][j];
}

void add9x9(double A[9][9], double B[9][9], double C[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            C[i][j] = A[i][j] + B[i][j];
}

void mult9x9_9x1(double A[9][9], double b[9], double c[9]) {
    for (int i = 0; i < 9; ++i) {
        c[i] = 0.0;
        for (int j = 0; j < 9; ++j)
            c[i] += A[i][j] * b[j];
    }
}

void transpose9x9(double A[9][9], double AT[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            AT[i][j] = A[j][i];
}

// ============ 4-BAR MECHANICS ============
void fourbarC(const double q[9], const SystemParams& p, double C[8]) {
    double Rx1 = q[0], Ry1 = q[1], th1 = q[2];
    double Rx2 = q[3], Ry2 = q[4], th2 = q[5];
    double Rx3 = q[6], Ry3 = q[7], th3 = q[8];
    
    double L1_2 = p.L1 / 2.0;
    double L2_2 = p.L2 / 2.0;
    double L3_2 = p.L3 / 2.0;
    
    double s1 = std::sin(th1), c1 = std::cos(th1);
    double s2 = std::sin(th2), c2 = std::cos(th2);
    double s3 = std::sin(th3), c3 = std::cos(th3);
    
    C[0] = Rx1 - L1_2 * c1;
    C[1] = Ry1 - L1_2 * s1;
    C[2] = Rx1 + L1_2 * s1 - Rx2 - L2_2 * s2;
    C[3] = Ry1 - L1_2 * c1 - Ry2 + L2_2 * c2;
    C[4] = Rx2 - L2_2 * s2 - Rx3 + L3_2 * s3;
    C[5] = Ry2 + L2_2 * c2 - Ry3 - L3_2 * c3;
    C[6] = Rx3 + L3_2 * s3;
    C[7] = Ry3 - L3_2 * c3;
}

void fourbarCq(const double q[9], const SystemParams& p, double Cq[8][9]) {
    double th1 = q[2], th2 = q[5], th3 = q[8];
    
    double L1_2 = p.L1 / 2.0;
    double L2_2 = p.L2 / 2.0;
    double L3_2 = p.L3 / 2.0;
    
    double s1 = std::sin(th1), c1 = std::cos(th1);
    double s2 = std::sin(th2), c2 = std::cos(th2);
    double s3 = std::sin(th3), c3 = std::cos(th3);
    
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 9; ++j)
            Cq[i][j] = 0.0;
    
    Cq[0][0] = 1.0;          Cq[0][2] = L1_2 * s1;
    Cq[1][1] = 1.0;          Cq[1][2] = -(L1_2 * c1);
    Cq[2][0] = 1.0;          Cq[2][2] = -(L1_2 * s1);   Cq[2][3] = -1.0;  Cq[2][5] = -(L2_2 * c2);
    Cq[3][1] = 1.0;          Cq[3][2] = L1_2 * c1;      Cq[3][4] = -1.0;  Cq[3][5] = L2_2 * c2;
    Cq[4][3] = 1.0;          Cq[4][5] = -(L2_2 * c2);   Cq[4][6] = -1.0;  Cq[4][8] = L3_2 * c3;
    Cq[5][4] = 1.0;          Cq[5][5] = L2_2 * c2;      Cq[5][7] = -1.0;  Cq[5][8] = -(L3_2 * c3);
    Cq[6][6] = 1.0;          Cq[6][8] = L3_2 * s3;
    Cq[7][7] = 1.0;          Cq[7][8] = -(L3_2 * c3);
}

void fourbarM(const SystemParams& p, double M[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            M[i][j] = 0.0;
    
    M[0][0] = p.m1;  M[1][1] = p.m1;  M[2][2] = p.I1;
    M[3][3] = p.m2;  M[4][4] = p.m2;  M[5][5] = p.I2;
    M[6][6] = p.m3;  M[7][7] = p.m3;  M[8][8] = p.I3;
}

void fourbarQe(double t, const SystemParams& p, double Qe[9]) {
    for (int i = 0; i < 9; ++i)
        Qe[i] = 0.0;
    
    Qe[1] = -p.m1 * p.g;
    Qe[4] = -p.m2 * p.g;
    Qe[7] = -p.m3 * p.g;
    Qe[2] = p.T0 * std::sin(p.win * t);
}

bool gaussianElimination(double A[17][17], double b[17], double sol[17]) {
    double a[17][17];
    double rhs[17];
    
    for (int i = 0; i < 17; ++i) {
        for (int j = 0; j < 17; ++j)
            a[i][j] = A[i][j];
        rhs[i] = b[i];
    }
    
    for (int col = 0; col < 17; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 17; ++row)
            if (std::fabs(a[row][col]) > std::fabs(a[pivot][col]))
                pivot = row;
        
        if (pivot != col) {
            for (int j = 0; j < 17; ++j)
                std::swap(a[col][j], a[pivot][j]);
            std::swap(rhs[col], rhs[pivot]);
        }
        
        if (std::fabs(a[col][col]) < 1e-12) {
            return false;
        }
        
        for (int row = col + 1; row < 17; ++row) {
            double factor = a[row][col] / a[col][col];
            for (int j = col; j < 17; ++j)
                a[row][j] -= factor * a[col][j];
            rhs[row] -= factor * rhs[col];
        }
    }
    
    for (int i = 16; i >= 0; --i) {
        sol[i] = rhs[i];
        for (int j = i + 1; j < 17; ++j)
            sol[i] -= a[i][j] * sol[j];
        sol[i] /= a[i][i];
    }
    
    return true;
}

void computeAccelerations(const State& state, const SystemParams& p, double accel[9]) {
    double q[9], dq[9];
    for (int i = 0; i < 9; ++i) {
        q[i] = state.q[i];
        dq[i] = state.dq[i];
    }
    
    double M[9][9];
    fourbarM(p, M);
    
    double Cq[8][9];
    fourbarCq(q, p, Cq);
    
    double Qe[9];
    fourbarQe(0.0, p, Qe);
    double Qv[9] = {0};
    
    double eps = 1e-6;
    double dCdqdq[8][9];
    
    for (int j = 0; j < 9; ++j) {
        double qp[9], qm[9];
        for (int k = 0; k < 9; ++k) {
            qp[k] = q[k];
            qm[k] = q[k];
        }
        qp[j] += eps;
        qm[j] -= eps;
        
        double Cqp[8][9], Cqm[8][9];
        fourbarCq(qp, p, Cqp);
        fourbarCq(qm, p, Cqm);
        
        for (int i = 0; i < 8; ++i) {
            double Cp_qdot = 0, Cm_qdot = 0;
            for (int k = 0; k < 9; ++k) {
                Cp_qdot += Cqp[i][k] * dq[k];
                Cm_qdot += Cqm[i][k] * dq[k];
            }
            dCdqdq[i][j] = (Cp_qdot - Cm_qdot) / (2.0 * eps);
        }
    }
    
    double gamma[8];
    for (int i = 0; i < 8; ++i) {
        gamma[i] = 0.0;
        for (int j = 0; j < 9; ++j)
            gamma[i] -= dCdqdq[i][j] * dq[j];
    }
    
    double A[17][17];
    double b[17];
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            A[i][j] = M[i][j];
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 8; ++j)
            A[i][9 + j] = Cq[j][i];
    
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 9; ++j)
            A[9 + i][j] = Cq[i][j];
    
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            A[9 + i][9 + j] = 0.0;
    
    for (int i = 0; i < 9; ++i)
        b[i] = Qe[i] + Qv[i];
    for (int i = 0; i < 8; ++i)
        b[9 + i] = gamma[i];
    
    double sol[17];
    if (!gaussianElimination(A, b, sol)) {
        return;
    }
    
    for (int i = 0; i < 9; ++i)
        accel[i] = sol[i];
}

void rk4Step(State& state, double dt, const SystemParams& p) {
    State k[4];
    double acc[9];
    
    computeAccelerations(state, p, acc);
    for (int i = 0; i < 9; ++i) {
        k[0].q[i]  = state.dq[i];
        k[0].dq[i] = acc[i];
    }
    
    State temp = state;
    for (int i = 0; i < 9; ++i) {
        temp.q[i]  += 0.5 * dt * k[0].q[i];
        temp.dq[i] += 0.5 * dt * k[0].dq[i];
    }
    computeAccelerations(temp, p, acc);
    for (int i = 0; i < 9; ++i) {
        k[1].q[i]  = temp.dq[i];
        k[1].dq[i] = acc[i];
    }
    
    temp = state;
    for (int i = 0; i < 9; ++i) {
        temp.q[i]  += 0.5 * dt * k[1].q[i];
        temp.dq[i] += 0.5 * dt * k[1].dq[i];
    }
    computeAccelerations(temp, p, acc);
    for (int i = 0; i < 9; ++i) {
        k[2].q[i]  = temp.dq[i];
        k[2].dq[i] = acc[i];
    }
    
    temp = state;
    for (int i = 0; i < 9; ++i) {
        temp.q[i]  += dt * k[2].q[i];
        temp.dq[i] += dt * k[2].dq[i];
    }
    computeAccelerations(temp, p, acc);
    for (int i = 0; i < 9; ++i) {
        k[3].q[i]  = temp.dq[i];
        k[3].dq[i] = acc[i];
    }
    
    for (int i = 0; i < 9; ++i) {
        state.q[i]  += (dt / 6.0) * (k[0].q[i] + 2*k[1].q[i] + 2*k[2].q[i] + k[3].q[i]);
        state.dq[i] += (dt / 6.0) * (k[0].dq[i] + 2*k[1].dq[i] + 2*k[2].dq[i] + k[3].dq[i]);
    }
}

// ============ KINEMATIC CLOSURE ============
void closureFromTheta1(double th1, const SystemParams& p,
                       double& Rx1, double& Ry1, double& th2,
                       double& Rx2, double& Ry2,
                       double& Rx3, double& Ry3, double& th3) {
    
    double Bx = p.L1 * std::cos(th1);
    double By = p.L1 * std::sin(th1);
    double Dx = p.L4, Dy = 0.0;
    
    double d = std::sqrt((Dx - Bx) * (Dx - Bx) + (Dy - By) * (Dy - By));
    
    if (d > (p.L2 + p.L3) || d < std::fabs(p.L2 - p.L3)) {
        return;
    }
    
    double ex_x = (Dx - Bx) / d;
    double ex_y = (Dy - By) / d;
    double ey_x = -ex_y, ey_y = ex_x;
    
    double a = (p.L2 * p.L2 - p.L3 * p.L3 + d * d) / (2.0 * d);
    double h_sq = p.L2 * p.L2 - a * a;
    if (h_sq < 0) h_sq = 0;
    double h = std::sqrt(h_sq);
    
    double P2_x = Bx + a * ex_x;
    double P2_y = By + a * ex_y;
    
    double C1_x = P2_x + h * ey_x;
    double C1_y = P2_y + h * ey_y;
    double C2_x = P2_x - h * ey_x;
    double C2_y = P2_y - h * ey_y;
    
    double Cx, Cy;
    if (C1_y >= C2_y) {
        Cx = C1_x; Cy = C1_y;
    } else {
        Cx = C2_x; Cy = C2_y;
    }
    
    th2 = std::atan2(Cy - By, Cx - Bx);
    th3 = std::atan2(Cy - Dy, Cx - Dx);
    
    Rx1 = (0.0 + Bx) / 2.0;
    Ry1 = (0.0 + By) / 2.0;
    Rx2 = (Bx + Cx) / 2.0;
    Ry2 = (By + Cy) / 2.0;
    Rx3 = (Cx + Dx) / 2.0;
    Ry3 = (Cy + Dy) / 2.0;
}

// ============ CORRECTED INDIRECT KALMAN FILTER ============

void indirectEKFInitialize(IndirectEKFState& ekf) {
    for (int i = 0; i < 9; ++i)
        ekf.delta_x[i] = 0.0;
    
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j)
            ekf.P[i][j] = 0.0;
    }
    
    for (int i = 0; i < 3; ++i)
        ekf.P[i][i] = 1e-6;
    
    for (int i = 3; i < 6; ++i)
        ekf.P[i][i] = 1e-6;
    
    for (int i = 6; i < 9; ++i)
        ekf.P[i][i] = 1e-8;
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            ekf.Q[i][j] = 0.0;
    
    // ========== TUNING: PRIORITIZE MATLAB OVER C++ ==========
    
    double sigma_accel_noise = 1.0;
    for (int i = 6; i < 9; ++i)
        ekf.Q[i][i] = sigma_accel_noise * sigma_accel_noise;
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            ekf.R[i][j] = 0.0;
    
    
    double sigma_angle = 0.00005;   // Very small: by assuming MATLAB angle is accurate
    double sigma_omega = 0.0005;    // Very small: by assuming MATLAB velocity is accurate
    
    ekf.R[2][2] = sigma_angle * sigma_angle;
    ekf.R[5][5] = sigma_omega * sigma_omega;
}

void computeNumericalJacobianStateTransition(
    const State& state, const SystemParams& p, double dt,
    double Phi[9][9]) {
    
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            Phi[i][j] = 0.0;
    
    for (int i = 0; i < 9; ++i)
        Phi[i][i] = 1.0;
    
    for (int i = 0; i < 3; ++i)
        Phi[i][i + 3] = dt;
    
    for (int i = 3; i < 6; ++i)
        Phi[i][i + 3] = dt;
}

void indirectEKFPredict(IndirectEKFState& ekf, const State& state, 
                        const SystemParams& p, double dt) {
    
    double Phi[9][9];
    computeNumericalJacobianStateTransition(state, p, dt, Phi);
    
    double Phi_P[9][9], Phi_P_PhiT[9][9];
    mult9x9(Phi, ekf.P, Phi_P);
    
    double Phi_T[9][9];
    transpose9x9(Phi, Phi_T);
    mult9x9(Phi_P, Phi_T, Phi_P_PhiT);
    
    add9x9(Phi_P_PhiT, ekf.Q, ekf.P);
    
    for (int i = 0; i < 9; ++i)
        ekf.delta_x[i] = 0.0;
}

void indirectEKFUpdate(IndirectEKFState& ekf, const State& state_nominal,
                       const State& state_meas, const SystemParams& p) {
    
    // ========== 2D MEASUREMENT MODEL ==========
    double z_meas[2] = {state_meas.q[2], state_meas.dq[2]};
    double z_pred[2] = {state_nominal.q[2], state_nominal.dq[2]};
    
    double y[2];
    y[0] = z_meas[0] - z_pred[0];
    y[1] = z_meas[1] - z_pred[1];
    
    // ========== INNOVATION COVARIANCE S (2×2) ========== need for double check 
    double S[2][2];
    S[0][0] = ekf.P[2][2] + ekf.R[2][2];
    S[0][1] = ekf.P[2][5];
    S[1][0] = ekf.P[5][2];
    S[1][1] = ekf.P[5][5] + ekf.R[5][5];
    
    // ========== INVERT S (2×2) ANALYTICALLY ========== need for double check 
    double det_S = S[0][0] * S[1][1] - S[0][1] * S[1][0];
    if (std::fabs(det_S) < 1e-14) {
        return;
    }
    
    double S_inv[2][2];
    S_inv[0][0] =  S[1][1] / det_S;
    S_inv[0][1] = -S[0][1] / det_S;
    S_inv[1][0] = -S[1][0] / det_S;
    S_inv[1][1] =  S[0][0] / det_S;
    
    // ========== KALMAN GAIN K (9×2) ========== need for doucle check
    double K[9][2];
    for (int i = 0; i < 9; ++i) {
        K[i][0] = ekf.P[i][2] * S_inv[0][0] + ekf.P[i][5] * S_inv[1][0];
        K[i][1] = ekf.P[i][2] * S_inv[0][1] + ekf.P[i][5] * S_inv[1][1];
    }
    
    // ========== ERROR STATE UPDATE ==========
    for (int i = 0; i < 9; ++i) {
        ekf.delta_x[i] = K[i][0] * y[0] + K[i][1] * y[1];
    }
    
    // ========== COVARIANCE UPDATE (Joseph Form) ==========
    double I_KH[9][9];
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            I_KH[i][j] = (i == j ? 1.0 : 0.0);
            I_KH[i][j] -= K[i][0] * (j == 2 ? 1.0 : 0.0);
            I_KH[i][j] -= K[i][1] * (j == 5 ? 1.0 : 0.0);
        }
    }
    
    double P_temp[9][9];
    mult9x9(I_KH, ekf.P, P_temp);
    
    double I_KH_T[9][9];
    transpose9x9(I_KH, I_KH_T);
    
    double P_joseph[9][9];
    mult9x9(P_temp, I_KH_T, P_joseph);
    
    double K_R[9][2];
    for (int i = 0; i < 9; ++i) {
        K_R[i][0] = K[i][0] * ekf.R[2][2];
        K_R[i][1] = K[i][1] * ekf.R[5][5];
    }
    
    double K_R_KT[9][9];
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            K_R_KT[i][j] = K_R[i][0] * K[j][0] + K_R[i][1] * K[j][1];
            ekf.P[i][j] = P_joseph[i][j] + K_R_KT[i][j];
        }
    }
}

void applyErrorCorrection(State& state, const IndirectEKFState& ekf,
                          const SystemParams& p) {
    
    // ========== STEP 1: CORRECT DRIVING LINK POSITION & VELOCITY ==========
    double theta1_corrected = state.q[2] + ekf.delta_x[2];
    double omega1_corrected = state.dq[2] + ekf.delta_x[5];
    
    // ========== STEP 2: RE-CLOSE KINEMATIC LOOP (POSITION) ==========
    double Rx1, Ry1, th2, Rx2, Ry2, Rx3, Ry3, th3;
    closureFromTheta1(theta1_corrected, p, 
                      Rx1, Ry1, th2, Rx2, Ry2, Rx3, Ry3, th3);
    
    state.q[0] = Rx1;
    state.q[1] = Ry1;
    state.q[2] = theta1_corrected;
    state.q[3] = Rx2;
    state.q[4] = Ry2;
    state.q[5] = th2;
    state.q[6] = Rx3;
    state.q[7] = Ry3;
    state.q[8] = th3;
    
    // ========== STEP 3: APPLY CORRECTIONS TO ALL VELOCITY STATES ==========
    state.dq[2] = omega1_corrected;
    state.dq[5] += ekf.delta_x[8];
    state.dq[8] += ekf.delta_x[8];
    
    // ========== STEP 4: RE-SOLVE CONSTRAINT VELOCITIES ==========
    double s1 = std::sin(theta1_corrected), c1 = std::cos(theta1_corrected);
    double s2 = std::sin(th2), c2 = std::cos(th2);
    double s3 = std::sin(th3), c3 = std::cos(th3);
    
    double a11 = p.L2 * c2,   a12 = -p.L3 * c3;
    double a21 = p.L2 * s2,   a22 = -p.L3 * s3;
    double det = a11 * a22 - a12 * a21;
    
    if (std::fabs(det) > 1e-8) {
        double rhs1 = -p.L1 * omega1_corrected * c1;
        double rhs2 = -p.L1 * omega1_corrected * s1;
        
        double omega2_corrected = (rhs1 * a22 - rhs2 * a12) / det;
        double omega3_corrected = (a11 * rhs2 - a21 * rhs1) / det;
        
        state.dq[5] = omega2_corrected;
        state.dq[8] = omega3_corrected;
    }
    
    // ========== STEP 5: UPDATE CARTESIAN VELOCITIES OF COM ==========
    double L1_2 = p.L1 / 2.0, L2_2 = p.L2 / 2.0, L3_2 = p.L3 / 2.0;
    
    state.dq[0] = -L1_2 * s1 * omega1_corrected;
    state.dq[1] =  L1_2 * c1 * omega1_corrected;
    
    state.dq[3] = state.dq[0] - L2_2 * s2 * state.dq[5];
    state.dq[4] = state.dq[1] + L2_2 * c2 * state.dq[5];
    
    state.dq[6] = -L3_2 * s3 * state.dq[8];
    state.dq[7] =  L3_2 * c3 * state.dq[8];
}

void saveToCSV() {
    std::ofstream csv("FourbarAL_IndirectIKF.csv");
    csv << "Time,ThetaMAT,OmegaMAT,ThetaAL_CPP,OmegaAL_CPP,ThetaCorrected,OmegaCorrected,CycleTime_us\n";
    csv << std::fixed << std::setprecision(8);
    
    for (size_t i = 0; i < logData.time.size(); ++i) {
        csv << logData.time[i] << ","
            << logData.theta_matlab[i] << ","
            << logData.omega_matlab[i] << ","
            << logData.theta_al[i] << ","
            << logData.omega_al[i] << ","
            << logData.theta_corrected[i] << ","
            << logData.omega_corrected[i] << ","
            << logData.cycle_time_us[i] << "\n";
    }
    csv.close();
    std::cout << "\nSaved " << logData.time.size() << " samples to FourbarAL_IndirectIKF.csv\n";
}

// ============ MAIN ============
int main() {
    std::cout << "\n" << std::string(140, '=') << "\n";
    std::cout << "C++ AUGMENTED LAGRANGE + INDIRECT KALMAN FILTER \n";
    std::cout << std::string(140, '=') << "\n\n";
    
    const double g = 9.81;
    const double L1 = 2.0, L2 = 8.0, L3 = 5.0, L4 = 10.2;
    const double m1 = 2.0, m2 = 8.0, m3 = 5.0;
    
    auto rodInertia = [](double m, double L) { return m * L * L / 12.0; };
    
    SystemParams p;
    p.g = g;
    p.L1 = L1; p.L2 = L2; p.L3 = L3; p.L4 = L4;
    p.m1 = m1; p.m2 = m2; p.m3 = m3;
    p.I1 = rodInertia(m1, L1);
    p.I2 = rodInertia(m2, L2);
    p.I3 = rodInertia(m3, L3);
    p.T0 = 0.0;
    p.win = 2.0 * PI * 0.02;
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Parameters: L=[" << L1 << ", " << L2 << ", " << L3 << ", " << L4 << "] m\n";
    std::cout << "            m=[" << m1 << ", " << m2 << ", " << m3 << "] kg\n";
    std::cout << "Integration: RK4 (dt=1ms)\n";
    std::cout << "  - Expected: Theta_corrected - Theta_MATLAB\n\n";
    
    const int serverPort = 5556;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
        std::cerr << "Socket creation failed\n";
        return 1;
    }
    
    int opt = 1;
    setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    
    struct sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(serverPort);
    addr.sin_addr.s_addr = INADDR_ANY;
    
    if (bind(sock, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        std::cerr << "Bind failed\n";
        return 1;
    }
    
    listen(sock, 1);
    std::cout << "Waiting for MATLAB client on port " << serverPort << "...\n";
    
    int client = accept(sock, NULL, NULL);
    if (client < 0) {
        std::cerr << "Accept failed\n";
        return 1;
    }
    std::cout << "MATLAB connected!\n\n";
    
    std::cout << std::string(180, '-') << "\n";
    std::cout << std::setw(8) << "Time(s)"
              << std::setw(14) << "Theta_MAT"
              << std::setw(14) << "Omega_MAT"
              << std::setw(14) << "Theta_AL"
              << std::setw(14) << "Omega_AL"
              << std::setw(14) << "Theta_EKF"
              << std::setw(14) << "Omega_EKF"
              << std::setw(14) << "Error(rad)"
              << std::setw(12) << "t_cycle(μs)\n";
    std::cout << std::string(180, '-') << "\n";
    
    double th1_0 = PI / 3.0;
    double Rx1_0, Ry1_0, th2_0, Rx2_0, Ry2_0, Rx3_0, Ry3_0, th3_0;
    closureFromTheta1(th1_0, p, Rx1_0, Ry1_0, th2_0, Rx2_0, Ry2_0, Rx3_0, Ry3_0, th3_0);
    
    State state;
    state.q[0] = Rx1_0;   state.q[1] = Ry1_0;   state.q[2] = th1_0;
    state.q[3] = Rx2_0;   state.q[4] = Ry2_0;   state.q[5] = th2_0;
    state.q[6] = Rx3_0;   state.q[7] = Ry3_0;   state.q[8] = th3_0;
    
    for (int i = 0; i < 9; ++i)
        state.dq[i] = 0.0;
    state.dq[2] = 0.5;
    
    IndirectEKFState ekf;
    indirectEKFInitialize(ekf);
    
    const double dt_tcp = 0.001;
    const int Nsamples = 10000;
    
    int samplecount = 0;
    double cycle_time_sum = 0.0;
    int num_cycles = 0;
    double max_cycle_time = 0.0;
    
    while (samplecount < Nsamples) {
        auto step_start = std::chrono::high_resolution_clock::now();
        
        double theta_meas, omega_meas;
        
        if (recv(client, &theta_meas, sizeof(double), 0) <= 0 ||
            recv(client, &omega_meas, sizeof(double), 0) <= 0) {
            std::cout << "\nConnection lost.\n";
            break;
        }
        
        theta_meas = swapDouble(theta_meas);
        omega_meas = swapDouble(omega_meas);
        
        samplecount++;
        double t = samplecount * dt_tcp;
        
        // ========== STEP 1: RK4 Integration ==========
        rk4Step(state, dt_tcp, p);
        
        double theta_al = state.q[2];
        double omega_al = state.dq[2];
        
        // ========== STEP 2: EKF Prediction ==========
        indirectEKFPredict(ekf, state, p, dt_tcp);
        
        // ========== STEP 3: EKF Update (MATLAB MEASUREMENT PULLS STATE) ==========
        State state_measured;
        for (int i = 0; i < 9; ++i) {
            state_measured.q[i] = state.q[i];
            state_measured.dq[i] = state.dq[i];
        }
        state_measured.q[2] = theta_meas;
        state_measured.dq[2] = omega_meas;
        
        indirectEKFUpdate(ekf, state, state_measured, p);
        
        // ========== STEP 4: Error Correction (WITH AGGRESSIVE KALMAN GAIN) ==========
        applyErrorCorrection(state, ekf, p);
        
        double theta_corrected = state.q[2];
        double omega_corrected = state.dq[2];
        
        // ========== Performance & Logging ==========
        auto step_end = std::chrono::high_resolution_clock::now();
        double cycle_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                   step_end - step_start).count();
        cycle_time_sum += cycle_time_us;
        num_cycles++;
        max_cycle_time = std::max(max_cycle_time, cycle_time_us);
        
        double error = theta_meas - theta_corrected;
        
        logData.time.push_back(t);
        logData.theta_matlab.push_back(theta_meas);
        logData.omega_matlab.push_back(omega_meas);
        logData.theta_al.push_back(theta_al);
        logData.omega_al.push_back(omega_al);
        logData.theta_corrected.push_back(theta_corrected);
        logData.omega_corrected.push_back(omega_corrected);
        logData.cycle_time_us.push_back(cycle_time_us);
        
        if (samplecount % 100 == 0) {
            printf("%8.3f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f %14.8f %12.2f\n",
                   t, theta_meas, omega_meas, theta_al, omega_al,
                   theta_corrected, omega_corrected, error, cycle_time_us);
        }
    }
    
    std::cout << std::string(180, '-') << "\n\n";
    
    close(client);
    close(sock);
    
    std::cout << "Simulation complete!\n";
    std::cout << "  Samples: " << samplecount << "\n";
    std::cout << "  Duration: " << (samplecount * dt_tcp) << " seconds\n";
    
    if (num_cycles > 0) {
        double avg_cycle_time = cycle_time_sum / num_cycles;
        std::cout << "\n  Performance:\n";
        std::cout << "    Average: " << std::fixed << std::setprecision(2)
                  << avg_cycle_time << " μs (" << (1e6 / avg_cycle_time) << " Hz)\n";
        std::cout << "    Maximum: " << max_cycle_time << " μs\n";
    }
    
    std::cout << "\n  Saving results...\n";
    saveToCSV();
    
    std::cout << "\nDone! Check FourbarAL_IndirectEKF.csv\n";    
    return 0;
}