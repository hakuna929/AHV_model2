clear; clc;

% 飞行器参数
p = ahv_params();

% =========================
% 1) 目标与初始位置（NED, m）
% =========================
H0 = 30e3;           % 初始高度 30 km
R0 = 8000e3;         % 初始距离 8000 km

% 飞行器初始坐标（NED）
x0 = 0; y0 = 0; z0 = -H0;

% 静止目标坐标（NED）——保证初始距离正好8000km
rt_n = [R0; 0; -H0];

% =========================
% 2) 初始姿态：对准目标LOS（零视线角误差）
% =========================
R_vec0 = rt_n - [x0; y0; z0];
dx = R_vec0(1); dy = R_vec0(2); dz = R_vec0(3);
Rhor = sqrt(dx^2 + dy^2);

psi0   = atan2(dy, dx);
theta0 = atan2(-dz, max(Rhor, 1)); % NED: z向下，所以向上用 -dz
phi0   = 0;

% =========================
% 3) 初始速度：大小由Ma0给，方向与机体x轴一致
%    因为机体已对准LOS，所以速度方向=LOS
% =========================
Ma0 = 6;
[a0, ~] = atmos_simple(H0);
V0 = Ma0 * a0;

u0 = V0; v0 = 0; w0 = 0;  % 关键：迎角/侧滑初值为0，航向由姿态决定

p0 = 0; q0 = 0; r0 = 0;
mf0 = p.mf;

X0 = [x0; y0; z0; u0; v0; w0; phi0; theta0; psi0; p0; q0; r0; mf0];

% =========================
% 4) 仿真设置：长航程 + 命中事件（50m）
% =========================
t_end_guess = R0 / max(V0,1);       % 粗略估计
tspan = [0, 1.5*t_end_guess];       % 给余量

opts = odeset('RelTol',1e-6,'AbsTol',1e-8, ...
              'Events', @(t,x) event_reach_target(t, x, rt_n, 50));

% =========================
% 5) 积分：每一步调用control_law得到U(t,X)
% =========================
[t, X] = ode45(@(t,x) ahv_6dof_ode(t, x, control_law(t, x, p, rt_n), p), ...
               tspan, X0, opts);

% =========================
% 6) 后处理：画距离等
% =========================
rm = X(:,1:3);
R  = rm - rt_n.';
range = sqrt(sum(R.^2,2));

V = sqrt(X(:,4).^2 + X(:,5).^2 + X(:,6).^2);

figure;
subplot(4,1,1);
plot(t, -X(:,3)/1000); grid on; ylabel('H (km)');

subplot(4,1,2);
plot(t, V); grid on; ylabel('V (m/s)');

subplot(4,1,3);
plot(t, X(:,8)*180/pi); grid on; ylabel('\theta (deg)');

subplot(4,1,4);
plot(t, range); grid on; ylabel('Range (m)'); xlabel('t (s)');
title('PN Guidance + Autopilot (stop when range < 50 m)');