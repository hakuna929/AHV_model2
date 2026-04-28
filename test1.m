clear; clc;

%% 时间
dt = 0.1;
T  = 3000;
t  = 0:dt:T;
N  = numel(t);

%% 地球/常量
Re = 6378137.0;
h0 = 30e3;
Rfix = Re + h0;                 % 固定半径球面

%% 目标（固定在同一球面）
latT = 13.35; lonT = 144.55;
rT = lla2ecef_cgcs2000(latT, lonT, h0); rT = rT(:);
rT = Rfix * rT / norm(rT);      % 强制到同一半径
vT = [0;0;0];

%% 初始
lat0 = 19.2; lon0 = 110.5;
r = lla2ecef_cgcs2000(lat0, lon0, h0); r = r(:);
r = Rfix * r / norm(r);

% 初速度：沿切平面朝目标投影
V0 = 1700;                       % m/s
u_up = r / norm(r);
dir0 = rT - (rT.'*u_up)*u_up;
dir0 = dir0 / norm(dir0);
v = V0 * dir0;

%% 导引参数（仅切向控制）
Npn = 2.5;
Kpp = 1.2;

a_tan_max = 120;                 % 切向机动能力
a_long_max = 20;                 % 速度调节
V_cmd_far = 1700;
V_cmd_mid = 1200;
V_cmd_near = 800;

%% 记录
Dist = zeros(N,1);
VcHist = zeros(N,1);
VHist = zeros(N,1);
Rhist = zeros(N,3);

hit_k = N;
for k = 1:N
    % 当前局部基
    u_up = r / norm(r);

    % 目标相对
    r_rel = rT - r;
    R = max(norm(r_rel),1);
    u_los = r_rel / R;

    % 速度
    V = max(norm(v),1e-6);
    u_v = v / V;

    % 目标相对速度
    v_rel = vT - v;
    Vc = -dot(r_rel, v_rel)/R;
    omega_los = cross(r_rel, v_rel)/(R^2 + 1e-9);

    %% 切向导引加速度（PP + PN）
    e_dir = u_los - u_v;
    a_pp = Kpp * V * e_dir;
    a_pn = Npn * max(Vc,0) * cross(u_los, omega_los);

    % 只保留切向分量（去径向）
    a_tan = a_pp + a_pn;
    a_tan = a_tan - dot(a_tan,u_up)*u_up;

    nt = norm(a_tan);
    if nt > a_tan_max
        a_tan = a_tan/nt * a_tan_max;
    end

    %% 速度调节（沿速度方向）
    if R > 800e3
        V_cmd = V_cmd_far;
    elseif R > 250e3
        V_cmd = V_cmd_mid;
    else
        V_cmd = V_cmd_near;
    end
    a_long = 0.03*(V_cmd - V);
    a_long = max(min(a_long, a_long_max), -a_long_max);
    a_long_vec = a_long * u_v;

    % 总控制（只用于方向/速度）
    a_cmd = a_tan + a_long_vec;

    %% 积分
    v = v + a_cmd*dt;
    r = r + v*dt;

    %% 关键：球面约束（防撞地）
    % 1) 位置投影回固定半径
    r = Rfix * r / norm(r);

    % 2) 速度去径向分量，保证始终切向飞行
    u_up = r / norm(r);
    v = v - dot(v,u_up)*u_up;

    %% 记录
    Dist(k) = R;
    VcHist(k)=Vc;
    VHist(k)=norm(v);
    Rhist(k,:)=r.';

    % 命中
    if R < 5e3
        hit_k = k;
        fprintf('HIT at t=%.1f s, miss=%.1f m\n', t(k), R);
        break;
    end

    if any(~isfinite([r;v;a_cmd]))
        hit_k = k;
        fprintf('NaN/Inf at step %d\n',k);
        break;
    end
end

idx = 1:hit_k;

%% 图1 距离
figure;
plot(t(idx), Dist(idx)/1000, 'm', 'LineWidth', 1.8); grid on;
xlabel('Time (s)'); ylabel('Distance (km)');
title('Spherical-Constraint Minimal Test: Distance');

%% 图2 闭合速度
figure;
plot(t(idx), VcHist(idx), 'b', 'LineWidth', 1.4); grid on; yline(0,'r--');
xlabel('Time (s)'); ylabel('V_c (m/s)');
title('Closing Speed');

%% 图3 速度
figure;
plot(t(idx), VHist(idx), 'k', 'LineWidth', 1.4); grid on;
xlabel('Time (s)'); ylabel('|v| (m/s)');
title('Speed');

%% 图4 球面轨迹
figure;
plot3(Rhist(idx,1)/1e3, Rhist(idx,2)/1e3, Rhist(idx,3)/1e3, 'b', 'LineWidth',1.4); hold on; grid on; axis equal;
[xe,ye,ze]=sphere(80);
surf(Rfix*xe/1e3, Rfix*ye/1e3, Rfix*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none');
plot3(rT(1)/1e3,rT(2)/1e3,rT(3)/1e3,'ro','MarkerFaceColor','r');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Trajectory on Fixed-Altitude Sphere');
legend('Vehicle','Altitude Sphere','Target');