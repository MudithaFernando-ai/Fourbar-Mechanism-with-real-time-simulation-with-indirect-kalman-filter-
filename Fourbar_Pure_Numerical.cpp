// Double checking is on going with thesis writing and reading
// Numerical Jacobian: Cq computed via finite differences on C(q)
// Analytical + Augmented Lagrange + IKF
// No torque is added
//
// Compile:
// g++ -std=c++11 -O3 Fourbar_numerical_10_2_2026.cpp -o fourbar_server
//
// Run:
// ./fourbar_server


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

// ====================== STRUCTURES ======================
struct SystemParams {
    double g, L1, L2, L3, L4;
    double m1, m2, m3;
    double I1, I2, I3;
    double T0, win;
};

struct State {
    double q[9];   // [Rx1, Ry1, θ1, Rx2, Ry2, θ2, Rx3, Ry3, θ3]
    double dq[9];  // velocities
};

struct IndirectEKFState {
    double delta_x[9];
    double P[9][9];
    double Q[9][9];
    double R[9][9];
};

struct LogData {
    std::vector<double> time;

    std::vector<double> theta_matlab, omega_matlab, alpha_matlab;
    std::vector<double> theta_cpp, omega_cpp, alpha_cpp;
    std::vector<double> theta_ikf, omega_ikf, alpha_ikf;
    std::vector<double> cycle_time_us;
};

LogData logData;

// ====================== BYTE SWAP ======================
double swapDouble(double value) {
    uint64_t temp;
    std::memcpy(&temp, &value, sizeof(double));
    temp = ((temp & 0xFF00000000000000ULL) >> 56) |
           ((temp & 0x00FF000000000000ULL) >> 40) |
           ((temp & 0x0000FF0000000000ULL) >> 24) |
           ((temp & 0x000000FF00000000ULL) >> 8)  |
           ((temp & 0x00000000FF000000ULL) << 8)  |
           ((temp & 0x0000000000FF0000ULL) << 24) |
           ((temp & 0x000000000000FF00ULL) << 40) |
           ((temp & 0x00000000000000FFULL) << 56);
    double result;
    std::memcpy(&result, &temp, sizeof(double));
    return result;
}

// ====================== TCP RECV ALL ======================
static ssize_t recvAll(int sockfd, void* buf, size_t len) {
    size_t total = 0;
    char* p = static_cast<char*>(buf);
    while (total < len) {
        ssize_t n = recv(sockfd, p + total, len - total, 0);
        if (n <= 0) return n;
        total += static_cast<size_t>(n);
    }
    return static_cast<ssize_t>(total);
}

// ====================== MATRIX UTILITIES ======================
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

void transpose9x9(double A[9][9], double AT[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            AT[i][j] = A[j][i];
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
        if (std::fabs(a[col][col]) < 1e-12)
            return false;
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

// ====================== 4-BAR CONSTRAINTS (NEW) ======================
// Constraint vector C(q)=0 (8 constraints, 9 coordinates)
void fourbarC(const double q[9], const SystemParams& p, double C[8]) {
    const double Rx1 = q[0], Ry1 = q[1], th1 = q[2];
    const double Rx2 = q[3], Ry2 = q[4], th2 = q[5];
    const double Rx3 = q[6], Ry3 = q[7], th3 = q[8];

    const double L1_2 = p.L1 / 2.0;
    const double L2_2 = p.L2 / 2.0;
    const double L3_2 = p.L3 / 2.0;

    // Ground pin at A = (0,0) to link1 COM
    C[0] = Rx1 - L1_2 * std::cos(th1);
    C[1] = Ry1 - L1_2 * std::sin(th1);

    // Joint B: end of link1 equals left end of link2
    C[2] = (Rx1 + L1_2 * std::cos(th1)) - (Rx2 - L2_2 * std::cos(th2));
    C[3] = (Ry1 + L1_2 * std::sin(th1)) - (Ry2 - L2_2 * std::sin(th2));

    // Joint C: right end of link2 equals left end of link3
    // (Signs chosen to match your original analytical Cq pattern)
    C[4] = (Rx2 + L2_2 * std::cos(th2)) - (Rx3 + L3_2 * std::cos(th3));
    C[5] = (Ry2 + L2_2 * std::sin(th2)) - (Ry3 + L3_2 * std::sin(th3));

    // Ground pin at D = (L4,0) to link3 COM
    C[6] = Rx3 - (p.L4 + L3_2 * std::cos(th3));
    C[7] = Ry3 - (0.0  + L3_2 * std::sin(th3));
}

// Numerical Jacobian Cq = dC/dq via central differences
void fourbarCq_numerical(const double q[9], const SystemParams& p, double Cq[8][9]) {
    const double eps = 1e-6;

    // Base constraints
    double C0[8];
    fourbarC(q, p, C0);

    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 9; ++j)
            Cq[i][j] = 0.0;

    for (int j = 0; j < 9; ++j) {
        double qp[9], qm[9];
        for (int k = 0; k < 9; ++k) {
            qp[k] = q[k];
            qm[k] = q[k];
        }

        double h = eps * (std::fabs(q[j]) + 1.0);
        qp[j] += h;
        qm[j] -= h;

        double Cp[8], Cm[8];
        fourbarC(qp, p, Cp);
        fourbarC(qm, p, Cm);

        for (int i = 0; i < 8; ++i) {
            Cq[i][j] = (Cp[i] - Cm[i]) / (2.0 * h);
        }
    }
}

// ====================== MASS + FORCES ======================
void fourbarM(const SystemParams& p, double M[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            M[i][j] = 0.0;
    M[0][0] = p.m1; M[1][1] = p.m1; M[2][2] = p.I1;
    M[3][3] = p.m2; M[4][4] = p.m2; M[5][5] = p.I2;
    M[6][6] = p.m3; M[7][7] = p.m3; M[8][8] = p.I3;
}

void fourbarQe(double t, const SystemParams& p, double Qe[9]) {
    for (int i = 0; i < 9; ++i)
        Qe[i] = 0.0;
    Qe[1] = -p.m1 * p.g;
    Qe[4] = -p.m2 * p.g;
    Qe[7] = -p.m3 * p.g;
    Qe[2] = p.T0 * std::sin(p.win * t);
}

// ====================== AUGMENTED/LAGRANGE SOLVE ======================
void computeAccelerations(const State& state, const SystemParams& p, double accel[9]) {
    double q[9], dq[9];
    for (int i = 0; i < 9; ++i) {
        q[i]  = state.q[i];
        dq[i] = state.dq[i];
    }

    double M[9][9];
    fourbarM(p, M);

    // (CHANGED) numerical constraint Jacobian
    double Cq[8][9];
    fourbarCq_numerical(q, p, Cq);

    double Qe[9];
    fourbarQe(0.0, p, Qe);
    double Qv[9] = {0};

    // gamma = -d/dt(Cq*dq) approximated numerically as in your original code
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
        fourbarCq_numerical(qp, p, Cqp);
        fourbarCq_numerical(qm, p, Cqm);

        for (int i = 0; i < 8; ++i) {
            double Cp_qdot = 0.0, Cm_qdot = 0.0;
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
    if (!gaussianElimination(A, b, sol))
        return;

    for (int i = 0; i < 9; ++i)
        accel[i] = sol[i];
}

// ====================== RK4 ======================
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
        state.q[i]  += (dt / 6.0) * (k[0].q[i]  + 2*k[1].q[i]  + 2*k[2].q[i]  + k[3].q[i]);
        state.dq[i] += (dt / 6.0) * (k[0].dq[i] + 2*k[1].dq[i] + 2*k[2].dq[i] + k[3].dq[i]);
    }
}

// ====================== CLOSURE FROM θ1 ======================
void closureFromTheta1(double th1, const SystemParams& p,
                       double& Rx1, double& Ry1, double& th2,
                       double& Rx2, double& Ry2,
                       double& Rx3, double& Ry3, double& th3) {
    double Ax = 0.0, Ay = 0.0;
    double Dx = p.L4, Dy = 0.0;
    double Bx = Ax + p.L1 * std::cos(th1);
    double By = Ay + p.L1 * std::sin(th1);

    double dx = Dx - Bx;
    double dy = Dy - By;
    double d  = std::sqrt(dx*dx + dy*dy);

    if (d > (p.L2 + p.L3) || d < std::fabs(p.L2 - p.L3)) {
        th2 = th1;
        th3 = th1;
        Rx1 = (p.L1/2.0)*std::cos(th1);
        Ry1 = (p.L1/2.0)*std::sin(th1);
        Rx2 = Rx1;
        Ry2 = Ry1;
        Rx3 = p.L4 + (p.L3/2.0)*std::cos(th1);
        Ry3 = (p.L3/2.0)*std::sin(th1);
        return;
    }

    double ex_x = dx / d;
    double ex_y = dy / d;
    double ey_x = -ex_y;
    double ey_y =  ex_x;

    double a = (p.L2*p.L2 - p.L3*p.L3 + d*d)/(2.0*d);
    double h_sq = p.L2*p.L2 - a*a;
    if (h_sq < 0.0) h_sq = 0.0;
    double h = std::sqrt(h_sq);

    double P2_x = Bx + a*ex_x;
    double P2_y = By + a*ex_y;

    double C1_x = P2_x + h*ey_x;
    double C1_y = P2_y + h*ey_y;
    double C2_x = P2_x - h*ey_x;
    double C2_y = P2_y - h*ey_y;

    double Cx, Cy;
    if (C1_y >= C2_y) { Cx = C1_x; Cy = C1_y; }
    else              { Cx = C2_x; Cy = C2_y; }

    th2 = std::atan2(Cy - By, Cx - Bx);
    th3 = std::atan2(Cy - Dy, Cx - Dx);

    Rx1 = (p.L1/2.0)*std::cos(th1);
    Ry1 = (p.L1/2.0)*std::sin(th1);

    Rx2 = 0.5*(Bx + Cx);
    Ry2 = 0.5*(By + Cy);

    Rx3 = p.L4 + (p.L3/2.0)*std::cos(th3);
    Ry3 =        (p.L3/2.0)*std::sin(th3);
}

// ====================== INDIRECT EKF ======================
void indirectEKFInitialize(IndirectEKFState& ekf) {
    for (int i = 0; i < 9; ++i)
        ekf.delta_x[i] = 0.0;

    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            ekf.P[i][j] = 0.0;

    for (int i = 0; i < 3; ++i)
        ekf.P[i][i] = 1e-6;
    for (int i = 3; i < 6; ++i)
        ekf.P[i][i] = 1e-6;
    for (int i = 6; i < 9; ++i)
        ekf.P[i][i] = 1e-8;

    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            ekf.Q[i][j] = 0.0;
    ekf.Q[2][2] = 1e-6;  // theta1
    ekf.Q[5][5] = 1e-4;  // omega1

    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            ekf.R[i][j] = 0.0;
    double sigma_angle = 0.001;
    double sigma_omega = 0.00001;
    ekf.R[2][2] = sigma_angle * sigma_angle;
    ekf.R[5][5] = sigma_omega * sigma_omega;
}

void computeNumericalJacobianStateTransition(
    const State& state, const SystemParams& p, double dt,
    double Phi[9][9]) {
    (void)state; (void)p;

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

void indirectEKFUpdate(IndirectEKFState& ekf, const State& est,
                       const State& meas, const SystemParams& p) {
    (void)p;
    double z_meas[2] = {meas.q[2], meas.dq[2]};
    double z_pred[2] = {est.q[2],  est.dq[2]};

    double y[2];
    y[0] = z_meas[0] - z_pred[0];
    y[1] = z_meas[1] - z_pred[1];

    double S[2][2];
    S[0][0] = ekf.P[2][2] + ekf.R[2][2];
    S[0][1] = ekf.P[2][5];
    S[1][0] = ekf.P[5][2];
    S[1][1] = ekf.P[5][5] + ekf.R[5][5];

    double det_S = S[0][0]*S[1][1] - S[0][1]*S[1][0];
    if (std::fabs(det_S) < 1e-14) return;

    double S_inv[2][2];
    S_inv[0][0] =  S[1][1] / det_S;
    S_inv[0][1] = -S[0][1] / det_S;
    S_inv[1][0] = -S[1][0] / det_S;
    S_inv[1][1] =  S[0][0] / det_S;

    double K[9][2];
    for (int i = 0; i < 9; ++i) {
        K[i][0] = ekf.P[i][2]*S_inv[0][0] + ekf.P[i][5]*S_inv[1][0];
        K[i][1] = ekf.P[i][2]*S_inv[0][1] + ekf.P[i][5]*S_inv[1][1];
    }

    for (int i = 0; i < 9; ++i)
        ekf.delta_x[i] = K[i][0]*y[0] + K[i][1]*y[1];

    double I_KH[9][9];
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) {
            I_KH[i][j] = (i == j ? 1.0 : 0.0);
            I_KH[i][j] -= K[i][0]*(j == 2 ? 1.0 : 0.0);
            I_KH[i][j] -= K[i][1]*(j == 5 ? 1.0 : 0.0);
        }

    double P_temp[9][9];
    mult9x9(I_KH, ekf.P, P_temp);

    double I_KH_T[9][9];
    transpose9x9(I_KH, I_KH_T);

    double P_joseph[9][9];
    mult9x9(P_temp, I_KH_T, P_joseph);

    double K_R[9][2];
    for (int i = 0; i < 9; ++i) {
        K_R[i][0] = K[i][0]*ekf.R[2][2];
        K_R[i][1] = K[i][1]*ekf.R[5][5];
    }

    double K_R_KT[9][9];
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) {
            K_R_KT[i][j] = K_R[i][0]*K[j][0] + K_R[i][1]*K[j][1];
            ekf.P[i][j]  = P_joseph[i][j] + K_R_KT[i][j];
        }
}

// ====================== APPLY ERROR CORRECTION ======================
void applyErrorCorrection(State& est, const IndirectEKFState& ekf,
                          const SystemParams& p) {
    double theta1_corrected = est.q[2] + ekf.delta_x[2];
    double omega1_corrected = est.dq[2] + ekf.delta_x[5];

    double Rx1, Ry1, th2, Rx2, Ry2, Rx3, Ry3, th3;
    closureFromTheta1(theta1_corrected, p,
                      Rx1, Ry1, th2, Rx2, Ry2, Rx3, Ry3, th3);

    est.q[0] = Rx1; est.q[1] = Ry1; est.q[2] = theta1_corrected;
    est.q[3] = Rx2; est.q[4] = Ry2; est.q[5] = th2;
    est.q[6] = Rx3; est.q[7] = Ry3; est.q[8] = th3;

    est.dq[2] = omega1_corrected;

    double s1 = std::sin(theta1_corrected), c1 = std::cos(theta1_corrected);
    double s2 = std::sin(th2),             c2 = std::cos(th2);
    double s3 = std::sin(th3),             c3 = std::cos(th3);

    double a11 = p.L2 * c2, a12 = -p.L3 * c3;
    double a21 = p.L2 * s2, a22 = -p.L3 * s3;
    double det = a11 * a22 - a12 * a21;

    if (std::fabs(det) > 1e-8) {
        double rhs1 = -p.L1 * omega1_corrected * c1;
        double rhs2 = -p.L1 * omega1_corrected * s1;

        double omega2_corrected = (rhs1 * a22 - rhs2 * a12) / det;
        double omega3_corrected = (a11 * rhs2 - a21 * rhs1) / det;

        est.dq[5] = omega2_corrected;
        est.dq[8] = omega3_corrected;
    }

    double L1_2 = p.L1 / 2.0, L2_2 = p.L2 / 2.0, L3_2 = p.L3 / 2.0;

    est.dq[0] = -L1_2 * s1 * omega1_corrected;
    est.dq[1] =  L1_2 * c1 * omega1_corrected;

    est.dq[3] = est.dq[0] - L2_2 * s2 * est.dq[5];
    est.dq[4] = est.dq[1] + L2_2 * c2 * est.dq[5];

    est.dq[6] = -L3_2 * s3 * est.dq[8];
    est.dq[7] =  L3_2 * c3 * est.dq[8];
}

// ====================== CSV SAVE ======================
void saveToCSV() {
    std::ofstream csv("FourbarAL_IKF_Numerical.csv");
    csv << "Time,"
           "ThetaMAT,OmegaMAT,AlphaMAT,"
           "ThetaCPP,OmegaCPP,AlphaCPP,"
           "ThetaIKF,OmegaIKF,AlphaIKF,"
           "CycleTime_us\n";
    csv << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < logData.time.size(); ++i) {
        csv << logData.time[i]          << ","
            << logData.theta_matlab[i]  << ","
            << logData.omega_matlab[i]  << ","
            << logData.alpha_matlab[i]  << ","
            << logData.theta_cpp[i]     << ","
            << logData.omega_cpp[i]     << ","
            << logData.alpha_cpp[i]     << ","
            << logData.theta_ikf[i]     << ","
            << logData.omega_ikf[i]     << ","
            << logData.alpha_ikf[i]     << ","
            << logData.cycle_time_us[i] << "\n";
    }
    csv.close();
    std::cout << "\nSaved " << logData.time.size()
              << " samples to FourbarAL_IKF_Numerical.csv\n";
}

// ====================== MAIN ======================
int main() {
    const double g  = 9.81;
    const double L1 = 2.0, L2 = 8.0, L3 = 5.0, L4 = 10.2;
    const double m1 = 2.0, m2 = 8.0, m3 = 5.0;

    auto rodInertia = [](double m, double L) { return m * L * L / 12.0; };

    SystemParams p;
    p.g  = g;
    p.L1 = L1; p.L2 = L2; p.L3 = L3; p.L4 = L4;
    p.m1 = m1; p.m2 = m2; p.m3 = m3;
    p.I1 = rodInertia(m1, L1);
    p.I2 = rodInertia(m2, L2);
    p.I3 = rodInertia(m3, L3);
    p.T0 = 0.0;
    p.win = 2.0 * PI * 0.02;

    const int serverPort = 5556;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) { std::cerr << "Socket creation failed\n"; return 1; }
    int opt = 1;
    setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    struct sockaddr_in addr;
    addr.sin_family      = AF_INET;
    addr.sin_port        = htons(serverPort);
    addr.sin_addr.s_addr = INADDR_ANY;
    if (bind(sock, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        std::cerr << "Bind failed\n"; return 1;
    }
    listen(sock, 1);
    std::cout << "Waiting for MATLAB client on port " << serverPort << "...\n";
    int client = accept(sock, NULL, NULL);
    if (client < 0) { std::cerr << "Accept failed\n"; return 1; }
    std::cout << "MATLAB connected!\n";

    double th1_0 = PI / 6.0;
    double Rx1_0, Ry1_0, th2_0, Rx2_0, Ry2_0, Rx3_0, Ry3_0, th3_0;
    closureFromTheta1(th1_0, p, Rx1_0, Ry1_0, th2_0, Rx2_0, Ry2_0, Rx3_0, Ry3_0, th3_0);

    State plant, est;
    plant.q[0] = Rx1_0; plant.q[1] = Ry1_0; plant.q[2] = th1_0;
    plant.q[3] = Rx2_0; plant.q[4] = Ry2_0; plant.q[5] = th2_0;
    plant.q[6] = Rx3_0; plant.q[7] = Ry3_0; plant.q[8] = th3_0;
    for (int i = 0; i < 9; ++i) plant.dq[i] = 0.0;

    est = plant;

    IndirectEKFState ekf;
    indirectEKFInitialize(ekf);

    const double dt_tcp   = 0.001;
    const int    Nsamples = 10000;

    int    samplecount    = 0;
    double max_cycle_time = 0.0;

    double omega_ikf_prev = 0.0;
    bool   first_ikf      = true;

    while (samplecount < Nsamples) {
        auto step_start = std::chrono::high_resolution_clock::now();

        // Receive theta, omega, alpha from MATLAB (3 doubles)
        double theta_meas, omega_meas, alpha_meas;
        if (recvAll(client, &theta_meas, sizeof(double)) <= 0 ||
            recvAll(client, &omega_meas, sizeof(double)) <= 0 ||
            recvAll(client, &alpha_meas, sizeof(double)) <= 0) {
            std::cout << "\nConnection lost.\n";
            break;
        }
        theta_meas = swapDouble(theta_meas);
        omega_meas = swapDouble(omega_meas);
        alpha_meas = swapDouble(alpha_meas);

        samplecount++;
        double t = samplecount * dt_tcp;

        // 1) Plant integration
        rk4Step(plant, dt_tcp, p);
        double theta_cpp = plant.q[2];
        double omega_cpp = plant.dq[2];
        double acc[9];
        computeAccelerations(plant, p, acc);
        double alpha_cpp = acc[2];

        // 2) MATLAB measurement (all from TCP/IP)
        double theta_matlab = theta_meas;
        double omega_matlab = omega_meas;
        double alpha_matlab = alpha_meas;

        // 3) EKF prediction on est
        indirectEKFPredict(ekf, est, p, dt_tcp);

        // 4) Measurement state for EKF
        State meas = est;
        meas.q[2]  = theta_matlab;
        meas.dq[2] = omega_matlab;

        // 5) EKF update
        indirectEKFUpdate(ekf, est, meas, p);

        // 6) Correction
        applyErrorCorrection(est, ekf, p);
        double theta_ikf = est.q[2];
        double omega_ikf = est.dq[2];
        double alpha_ikf = 0.0;
        if (!first_ikf) alpha_ikf = (omega_ikf - omega_ikf_prev) / dt_tcp;
        omega_ikf_prev = omega_ikf;
        first_ikf      = false;

        // 7) Perf + log
        auto   step_end      = std::chrono::high_resolution_clock::now();
        double cycle_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                   step_end - step_start).count();
        max_cycle_time = std::max(max_cycle_time, cycle_time_us);
        double error = theta_matlab - theta_ikf;

        logData.time.push_back(t);
        logData.theta_matlab.push_back(theta_matlab);
        logData.omega_matlab.push_back(omega_matlab);
        logData.alpha_matlab.push_back(alpha_matlab);
        logData.theta_cpp.push_back(theta_cpp);
        logData.omega_cpp.push_back(omega_cpp);
        logData.alpha_cpp.push_back(alpha_cpp);
        logData.theta_ikf.push_back(theta_ikf);
        logData.omega_ikf.push_back(omega_ikf);
        logData.alpha_ikf.push_back(alpha_ikf);
        logData.cycle_time_us.push_back(cycle_time_us);

        if (samplecount % 100 == 0) {
            printf("%8.3f %14.6f %14.6f %14.6f "
                   "%14.6f %14.6f %14.6f "
                   "%14.6f %14.6f %14.6f "
                   "%14.8f %12.2f\n",
                   t,
                   theta_matlab, omega_matlab, alpha_matlab,
                   theta_cpp,     omega_cpp,     alpha_cpp,
                   theta_ikf,     omega_ikf,     alpha_ikf,
                   error, cycle_time_us);
        }
    }

    close(client);
    close(sock);
    saveToCSV();
    return 0;
}
