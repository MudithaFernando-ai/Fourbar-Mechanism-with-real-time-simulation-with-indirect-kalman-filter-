% MATLAB fourbar mechanism simulation (WITH TCP/IP)
clear; clc; close all;

%% Parameters
g  = 8;
L1 = 2.0;  L2 = 8.0;  L3 = 5.0;  L4 = 10.2;
m1 = 2.0;  m2 = 8.0;  m3 = 5.0;
rodInertia = @(m,L) m*L^2/12;
params.g = g;
params.L1=L1; params.L2=L2; params.L3=L3; params.L4=L4;
params.m1=m1; params.m2=m2; params.m3=m3;
params.I1=rodInertia(m1,L1);
params.I2=rodInertia(m2,L2);
params.I3=rodInertia(m3,L3);
params.T0=0.0;
params.win=2*pi*0.02;

%% ── NOISE PARAMETERS  ────────────────────────────────────────────
rng(42);
sigma_theta = 1e-3;   % measurement noise std – theta  (rad)
sigma_omega = 5e-3;   % measurement noise std – omega  (rad/s)
sigma_alpha = 1e-2;   % measurement noise std – alpha  (rad/s²)

%% Time
dt   = 0.001;
tend = 10.0;
tsol = 0:dt:tend;
N    = numel(tsol);

%% TCP/IP CONFIGURATION
serverIP = "169.254.131.136";  %  C++ server IP
serverPort = 5556;             % server port

%% CONNECT TO C++ SERVER
try
    tcpipClient = tcpip(serverIP, serverPort, 'NetworkRole', 'client');
    tcpipClient.InputBufferSize = 8192;
    tcpipClient.OutputBufferSize = 8192;
    tcpipClient.Timeout = 30;
    fopen(tcpipClient);
    fprintf(' Connected to C++ Augmented Lagrange server!\n\n');
catch ME
    fprintf(' Connection failed: %s\n', ME.message);
    return;
end

%% Initial configuration (closure)
th1_0 = pi/2;
[B0, C0, th2_0, th3_0] = closure_from_th1(th1_0, params);
Rx1_0 = (L1/2)*cos(th1_0);
Ry1_0 = (L1/2)*sin(th1_0);
Rx2_0 = 0.5*(B0(1) + C0(1));
Ry2_0 = 0.5*(B0(2) + C0(2));
Rx3_0 = L4 + (L3/2)*cos(th3_0);
Ry3_0 =       (L3/2)*sin(th3_0);
q0  = [Rx1_0; Ry1_0; th1_0; Rx2_0; Ry2_0; th2_0; Rx3_0; Ry3_0; th3_0];
dq0 = zeros(9,1);
x0  = [q0; dq0];
xsol = zeros(N,18);
xsol(1,:) = x0';

%% --- Logging for "show" each 0.001 s transfer ---
theta_log = zeros(N,1);
omega_log = zeros(N,1);
alpha_log = zeros(N,1);

theta_meas_log = zeros(N,1);   % (Added)
omega_meas_log = zeros(N,1);   % (Added)
alpha_meas_log = zeros(N,1);   % (Added)

wall_log  = zeros(N,1);   % wall-clock time at each transfer
txdt_log  = zeros(N,1);   % wall-clock interval between transfers

% initial values (True)
theta_log(1) = x0(3);
omega_log(1) = x0(12);
dx0 = fourbarode_lagrange(tsol(1), x0, params);
alpha_log(1) = dx0(12);   % ddq(3)

% initial values (Noisy)
theta_meas_log(1) = theta_log(1) + sigma_theta*randn();
omega_meas_log(1) = omega_log(1) + sigma_omega*randn();
alpha_meas_log(1) = alpha_log(1) + sigma_alpha*randn();

wall_log(1)  = 0;
txdt_log(1)  = 0;
showEveryStep = true;     % prints every 0.001 s (will slow MATLAB a lot)
printStride   = 1;        % keep 1 to print every step, or set 10/100 to reduce

%% Main loop (RK4 integration + TCP transmission)
tic;                 % For wall-clock timing
prevWall = toc;
for k = 1:N-1
    t = tsol(k);
    x = xsol(k,:)';
    
    % Integrate one step RK4
    k1 = fourbarode_lagrange(t,        x,            params);
    k2 = fourbarode_lagrange(t+dt/2.0, x+dt/2.0*k1,  params);
    k3 = fourbarode_lagrange(t+dt/2.0, x+dt/2.0*k2,  params);
    k4 = fourbarode_lagrange(t+dt,     x+dt*k3,      params);
    x_next = x + dt/6.0*(k1 + 2*k2 + 2*k3 + k4);
    xsol(k+1,:) = x_next';
    
    % Extract TRUE theta, omega, alpha
    theta_true = x_next(3);    % q(3)
    omega_true = x_next(12);   % dq(3)
    dxdt_next = fourbarode_lagrange(t+dt, x_next, params);
    alpha_true = dxdt_next(12); % ddq(3) 

    % Add NOISE to create the measured signals to be sent
    theta_current = theta_true + sigma_theta*randn();
    omega_current = omega_true + sigma_omega*randn();
    alpha_current = alpha_true + sigma_alpha*randn();
    
    % Log both
    theta_log(k+1) = theta_true;
    omega_log(k+1) = omega_true;
    alpha_log(k+1) = alpha_true;

    theta_meas_log(k+1) = theta_current;
    omega_meas_log(k+1) = omega_current;
    alpha_meas_log(k+1) = alpha_current;
    
    nowWall = toc;
    wall_log(k+1) = nowWall;
    txdt_log(k+1) = nowWall - prevWall;
    prevWall = nowWall;
    
    % ===============================================
    % STEP 8: TCP TRANSMISSION (send noisy theta, omega, alpha)
    % ===============================================
    try
        % Send as 3 doubles in one write (24 bytes total)
        fwrite(tcpipClient, [theta_current; omega_current; alpha_current], 'double');
    catch
        fprintf('\n TCP transmission error\n');
        break;
    end
    
    % Optional console display every step
    if showEveryStep && (mod(k,printStride)==0)
        fprintf('t=%.3f  theta=%.6f  omega=%.6f  alpha=%.6f  tx_dt=%.6f s\n', ...
            t+dt, theta_current, omega_current, alpha_current, txdt_log(k+1));
    end
end

%% CLEANUP
try
    fclose(tcpipClient);
    delete(tcpipClient);
catch
end

%% Plot results (theta, omega, alpha vs simulation time)
figure('Name','Fourbar TX signals (True vs Noisy)');

subplot(3,1,1); 
plot(tsol, theta_log, 'b', 'LineWidth', 1.4); hold on;
plot(tsol, theta_meas_log, 'r', 'LineWidth', 0.8);
grid on; ylabel('\theta_1 (rad)'); title('\theta_1 (True vs Measured)');
legend('True', 'Measured', 'Location','best');

subplot(3,1,2); 
plot(tsol, omega_log, 'b', 'LineWidth', 1.4); hold on;
plot(tsol, omega_meas_log, 'r', 'LineWidth', 0.8);
grid on; ylabel('\omega_1 (rad/s)'); title('\omega_1 (True vs Measured)');
legend('True', 'Measured', 'Location','best');

subplot(3,1,3); 
plot(tsol, alpha_log, 'b', 'LineWidth', 1.4); hold on;
plot(tsol, alpha_meas_log, 'r', 'LineWidth', 0.8);
grid on; ylabel('\alpha_1 (rad/s^2)'); xlabel('t (s)'); title('\alpha_1 (True vs Measured)');
legend('True', 'Measured', 'Location','best');

figure('Name','Wall-clock transfer interval');
plot(tsol, txdt_log, 'LineWidth', 1.2); grid on;
xlabel('t (s)'); ylabel('Transfer interval (s)');

%% ===================== functions =====================
function dxdt = fourbarode_lagrange(t, x, p)
    q  = x(1:9);
    dq = x(10:18);
    M  = fourbarM(p);
    Cq = fourbarCq(q, p);
    Qe = fourbarQe(t, p);
    Qv = zeros(9,1);
    eps = 1e-6;
    nC  = size(Cq,1);
    nQ  = numel(q);
    dCdqdq = zeros(nC, nQ);
    for j = 1:nQ
        qp = q; qm = q;
        qp(j) = qp(j) + eps;
        qm(j) = qm(j) - eps;
        Cqp = fourbarCq(qp, p);
        Cqm = fourbarCq(qm, p);
        dCdqdq(:,j) = (Cqp*dq - Cqm*dq)/(2*eps);
    end
    gamma = -dCdqdq * dq;
    A = [M, Cq';
         Cq, zeros(nC,nC)];
    rhs = [Qe + Qv;
           gamma];
    sol = A \ rhs;
    ddq = sol(1:9);
    dxdt = zeros(18,1);
    dxdt(1:9)   = dq;
    dxdt(10:18) = ddq;
end

function M = fourbarM(p)
    M1 = diag([p.m1, p.m1, p.I1]);
    M2 = diag([p.m2, p.m2, p.I2]);
    M3 = diag([p.m3, p.m3, p.I3]);
    M  = blkdiag(M1, M2, M3);
end

function Qe = fourbarQe(t, p)
    Qe = zeros(9,1);
    Qe(2) = -p.m1*p.g;
    Qe(5) = -p.m2*p.g;
    Qe(8) = -p.m3*p.g;
    Qe(3) = Qe(3) + p.T0*sin(p.win*t);
    Qe(3) = Qe(3);
end

function Cq = fourbarCq(q, p)
    th1=q(3); th2=q(6); th3=q(9);
    L1=p.L1; L2=p.L2; L3=p.L3;
    s1=sin(th1); c1=cos(th1);
    s2=sin(th2); c2=cos(th2);
    s3=sin(th3); c3=cos(th3);
    Cq = zeros(8,9);
    Cq(1,1)=1;      Cq(1,3)= (L1/2)*s1;
    Cq(2,2)=1;      Cq(2,3)= -(L1/2)*c1;
    Cq(3,1)=1;      Cq(3,3)= -(L1/2)*s1;   Cq(3,4)=-1;  Cq(3,6)= -(L2/2)*s2;
    Cq(4,2)=1;      Cq(4,3)=  (L1/2)*c1;   Cq(4,5)=-1;  Cq(4,6)=  (L2/2)*c2;
    Cq(5,4)=1;      Cq(5,6)= -(L2/2)*s2;   Cq(5,7)=-1;  Cq(5,9)=  (L3/2)*s3;
    Cq(6,5)=1;      Cq(6,6)=  (L2/2)*c2;   Cq(6,8)=-1;  Cq(6,9)= -(L3/2)*c3;
    Cq(7,7)=1;      Cq(7,9)= (L3/2)*s3;
    Cq(8,8)=1;      Cq(8,9)= -(L3/2)*c3;
end

function [B, C, th2, th3] = closure_from_th1(th1, p)
    A = [0;0];
    D = [p.L4; 0];
    B = A + [p.L1*cos(th1); p.L1*sin(th1)];
    d = norm(D - B);
    if d > (p.L2 + p.L3) || d < abs(p.L2 - p.L3)
        error('No real closure');
    end
    ex = (D - B)/d;
    ey = [-ex(2); ex(1)];
    a = (p.L2^2 - p.L3^2 + d^2)/(2*d);
    h = sqrt(max(p.L2^2 - a^2, 0));
    P2 = B + a*ex;
    C1 = P2 + h*ey;
    C2 = P2 - h*ey;
    if C1(2) >= C2(2)
        C = C1;
    else
        C = C2;
    end
    th2 = atan2(C(2) - B(2), C(1) - B(1));
    th3 = atan2(C(2) - D(2), C(1) - D(1));
end

function [rA, rB, rC, rD] = joints_from_q(q, p)
    Rx1 = q(1); Ry1 = q(2); th1 = q(3);
    Rx2 = q(4); Ry2 = q(5); th2 = q(6);
    Rx3 = q(7); Ry3 = q(8); th3 = q(9);
    rA = [0; 0];
    rB = [Rx1 + (p.L1/2)*cos(th1); Ry1 + (p.L1/2)*sin(th1)];
    rC = [Rx2 + (p.L2/2)*cos(th2); Ry2 + (p.L2/2)*sin(th2)];
    rD = [p.L4; 0];
end
