%% Fourbar_9DOF_Animation.m
% MATLAB fourbar mechanism simulation (no TCP/IP, with plots + animation)
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

%% Time
dt   = 0.001;
tend = 10.0;
tsol = 0:dt:tend;
N    = numel(tsol);

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

%% Main loop (RK4 integration)
tic;
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
end
toc;

%% Plot theta1 and omega1
figure('Name', 'Four-bar Crank Kinematics');
subplot(2,1,1);
plot(tsol, xsol(:,3), 'b-', 'LineWidth', 2);
ylabel('$\theta_1$ (rad)');
title('Crank Angle vs Time');
grid on;

subplot(2,1,2);
plot(tsol, xsol(:,12), 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('$\omega_1$ (rad/s)');
title('Crank Angular Velocity vs Time');
grid on;

%% Animation (replay full simulation)
figure('Name', 'Four-bar Mechanism Animation');
fig = gcf;
fig.Position = [100, 100, 1000, 600];

% Precompute all joint positions for speed
rA_all = zeros(2,N); rB_all = zeros(2,N); rC_all = zeros(2,N); rD_all = zeros(2,N);
for k = 1:N
    [rA, rB, rC, rD] = joints_from_q(xsol(k,1:9)', params);
    rA_all(:,k) = rA; rB_all(:,k) = rB; rC_all(:,k) = rC; rD_all(:,k) = rD;
end

hAB = plot(NaN, NaN, 'b-', 'LineWidth', 4); hold on;
hBC = plot(NaN, NaN, 'g-', 'LineWidth', 4);
hCD = plot(NaN, NaN, 'r-', 'LineWidth', 4);
hGround = plot([0 L4], [0 0], 'k-', 'LineWidth', 6);
scatter([0 L4], [0 0], 200, 'ko', 'filled'); % Fixed joints A,D
hJoints = scatter(NaN, NaN, 150, 'k', 'filled'); % Moving joints B,C
hCoM1 = scatter(NaN, NaN, 100, 'b', 'filled'); % Crank CoM
hCoM2 = scatter(NaN, NaN, 100, 'g', 'filled'); % Coupler CoM
hCoM3 = scatter(NaN, NaN, 100, 'r', 'filled'); % Rocker CoM
title('9-DOF Four-bar Animation (L1=2, L2=8, L3=5, L4=10.2)');
xlabel('X (m)'); ylabel('Y (m)');
axis equal; grid on; xlim([-1 12]); ylim([-6 6]);
legend([hAB hBC hCD hGround], {'Crank AB', 'Coupler BC', 'Rocker CD', 'Ground'}, 'Location', 'best');

% Animate at 50 FPS (20ms/frame)
for k = 1:10:N  % Every 10th frame for smooth playback
    % Update links
    set(hAB, 'XData', [rA_all(1,k) rB_all(1,k)], 'YData', [rA_all(2,k) rB_all(2,k)]);
    set(hBC, 'XData', [rB_all(1,k) rC_all(1,k)], 'YData', [rB_all(2,k) rC_all(2,k)]);
    set(hCD, 'XData', [rC_all(1,k) rD_all(1,k)], 'YData', [rC_all(2,k) rD_all(2,k)]);
    
    % Update joints B,C
    set(hJoints, 'XData', [rB_all(1,k) rC_all(1,k)], 'YData', [rB_all(2,k) rC_all(2,k)]);
    
    % Update CoM positions
    com1 = [xsol(k,1); xsol(k,2)];
    com2 = [xsol(k,4); xsol(k,5)];
    com3 = [xsol(k,7); xsol(k,8)];
    set(hCoM1, 'XData', com1(1), 'YData', com1(2));
    set(hCoM2, 'XData', com2(1), 'YData', com2(2));
    set(hCoM3, 'XData', com3(1), 'YData', com3(2));
    
    % Status text
    text(0.02, 5.5, sprintf('t=%.2fs | θ₁=%.2f rad (%.1f°)', tsol(k), xsol(k,3), rad2deg(xsol(k,3))), ...
         'FontSize', 12, 'BackgroundColor', 'w');
    
    drawnow limitrate nocallbacks;
    pause(0.02);  % 50 FPS
end

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
    % Extract positions and angles from state q
    Rx1 = q(1); Ry1 = q(2); th1 = q(3);
    Rx2 = q(4); Ry2 = q(5); th2 = q(6);
    Rx3 = q(7); Ry3 = q(8); th3 = q(9);
    
    % Joint A (fixed origin)
    rA = [0; 0];
    
    % Joint B (end of crank AB)
    rB = [Rx1 + (p.L1/2)*cos(th1); Ry1 + (p.L1/2)*sin(th1)];
    
    % Joint C (end of coupler BC)
    rC = [Rx2 + (p.L2/2)*cos(th2); Ry2 + (p.L2/2)*sin(th2)];
    
    % Joint D (fixed ground end)
    rD = [p.L4; 0];
end
