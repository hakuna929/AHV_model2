% main_3dof_L1_terminal.m
% 3DOF终端版：L1航线跟踪 + 高度/速度保持 + 末段捕获增强
% 关键增强：
% 1) 末段自适应步长 dt_use（防跨步飞过）
% 2) 末段L1缩距 + 横向过载上调
% 3) 严格命中判据（R_hit）
% 4) 最近点判据（Vc<=0 且近距）防"飞过后还继续跑"
%
% 依赖函数：
% lla2ecef_cgcs2000, ecef2lla_cgcs2000, atmos_simple, aero_coeffs, thrust_coeffs

clear; clc;

%% ===================== 仿真时间 =====================
dt_base = 0.1;
T_end   = 3500;
N_max   = ceil(T_end/dt_base) * 5;  % 预留给末段小步长

%% ===================== 常量 =====================
mu  = 3.986004418e14;
Re  = 6378137.0;
we  = 7.2921150e-5;
omega_ie = [0;0;we];

use_simple_gravity = true;
use_rotation_terms = false;

%% ===================== 飞行器参数 =====================
m     = 671.33;
S_ref = 0.2986;

%% ===================== 初始/目标 =====================
h0 = 30e3; lat0 = 19.2; lon0 = 110.5;
r = lla2ecef_cgcs2000(lat0, lon0, h0); r = r(:);

latT = 13.35; lonT = 144.55; hT = 30e3;
rT = lla2ecef_cgcs2000(latT, lonT, hT); rT = rT(:);

[a0,~] = atmos_simple(h0);
V0 = 6.5 * a0;

% 航点（可扩展）
wps_lla = [lat0, lon0, h0;
           latT, lonT, hT];
nWP = size(wps_lla,1);
wps_ecef = zeros(nWP,3);
for i=1:nWP
    wps_ecef(i,:) = lla2ecef_cgcs2000(wps_lla(i,1), wps_lla(i,2), wps_lla(i,3));
end
wp_idx = 2;
r_wp_prev = wps_ecef(1,:)';
r_wp_next = wps_ecef(2,:)';
wp_switch_dist = 25e3;

% 初始速度沿航段切向（水平）
u_up0 = r/norm(r);
seg0 = r_wp_next - r_wp_prev;
seg0_h = seg0 - dot(seg0,u_up0)*u_up0;
if norm(seg0_h) < 1e-6, seg0_h = seg0; end
v = V0 * seg0_h / norm(seg0_h);

% 姿态（3DOF理想跟踪）
phi = 0; theta = 0; psi = 0;

%% ===================== 控制参数 =====================
% L1
L1_base  = 100e3;
L1_gainV = 25;
a_lat_max_cruise = 60;   % 巡航段
a_lat_max_term   = 140;  % 末段增强

% 纵向高度环
Kph = 0.010; Kih = 0.00002; Kdh = 0.08;
int_h = 0; int_h_lim = 8e4;
a_h_max = 25;

% alpha/phi 限幅
phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2*pi/180;
alpha_max = 18*pi/180;

% 速度控制（CT前馈+PI）
Kpv = 0.0018; Kiv = 0.00025;
int_v = 0; int_v_lim = 8e3;

CT_min = 0.00; CT_max = 0.45;
dT_min = 0.10; dT_max = 1.00;
tau_dT = 0.6; dT = 0.65;

% 马赫指令
M_cmd_far  = 6.2;
M_cmd_mid  = 5.2;
M_cmd_near = 4.5;
V_floor    = 750;

% 终端捕获参数
R_hit            = 1000;    % 1 km 命中
R_term1          = 120e3;   % 进入末段1
R_term2          = 40e3;    % 进入末段2（更激进）
R_pass_check     = 20e3;    % 最近点判据启用距离
R_abort_diverge  = 50e3;    % 若近距后发散可停

% 舵偏固定（3DOF简化）
de1 = 0; de2 = 0; dr = 0;

%% ===================== 记录数组 =====================
t_log = nan(N_max,1);
Rhist = nan(N_max,3);
Vhist = nan(N_max,3);
Hhist = nan(N_max,1);
Dist_hist = nan(N_max,1);
Vc_hist = nan(N_max,1);
Vmag_hist = nan(N_max,1);
Vcmd_hist = nan(N_max,1);
Mcmd_hist = nan(N_max,1);
dT_hist = nan(N_max,1);
CL_hist = nan(N_max,1);
CD_hist = nan(N_max,1);
CTff_hist = nan(N_max,1);
CTcmd_hist = nan(N_max,1);
L1_hist = nan(N_max,1);
alat_hist = nan(N_max,1);
wp_idx_hist = nan(N_max,1);

stop_reason = "completed";
k = 1;
t_now = 0;

prev_R = inf;
min_R = inf;

%% ===================== 主循环 =====================
while k <= N_max && t_now <= T_end
    % --- 当前地理量 ---
    [lat, lon, h] = ecef2lla_cgcs2000(r');
    if h <= 20e3
        stop_reason = "terrain impact";
        fprintf('Vehicle impacted terrain at t=%.1f s\n', t_now);
        break;
    end

    % --- 目标相对 ---
    r_rel_T = rT - r;
    R_to_T = norm(r_rel_T);
    u_los_T = r_rel_T / max(R_to_T,1);
    Vc_toward = dot(v, u_los_T);   % >0表示正在接近目标

    min_R = min(min_R, R_to_T);

    % --- 命中判据 ---
    if R_to_T <= R_hit
        stop_reason = "target reached";
        fprintf('Target reached at t=%.1f s, miss=%.1f m\n', t_now, R_to_T);
        break;
    end

    % --- 最近点判据：近距且闭合速度反号（飞过）---
    if (R_to_T < R_pass_check) && (Vc_toward <= 0)
        stop_reason = "passed closest approach";
        fprintf('Passed closest approach at t=%.1f s, min range=%.1f m\n', t_now, min_R);
        break;
    end

    % --- 近距后发散保护 ---
    if (prev_R < R_abort_diverge) && (R_to_T > prev_R + 30)  % 连续变大
        stop_reason = "diverging after near pass";
        fprintf('Diverging after near pass at t=%.1f s, range=%.1f m\n', t_now, R_to_T);
        break;
    end
    prev_R = R_to_T;

    % --- 自适应步长 ---
    if R_to_T < R_term2
        dt_use = 0.001;
    elseif R_to_T < R_term1
        dt_use = 0.01;
    else
        dt_use = dt_base;
    end

    % --- 大气 ---
    V = max(norm(v),1e-3);
    [a_snd, rho] = atmos_simple(max(h,0));
    Ma = max(min(V/max(a_snd,1e-3),8.0),0.0);
    qbar = 0.5*rho*V^2;

    % --- 局部基 ---
    u_up = r / norm(r);
    eE = [-sin(lon); cos(lon); 0];
    eN = [-sin(lat)*cos(lon); -sin(lat)*sin(lon); cos(lat)];

    v_h = v - dot(v,u_up)*u_up;
    Vh = max(norm(v_h),1e-6);
    u_vh = v_h / Vh;

    %% 航点切换
    if norm(r_wp_next - r) < wp_switch_dist && wp_idx < nWP
        wp_idx = wp_idx + 1;
        r_wp_prev = r_wp_next;
        r_wp_next = wps_ecef(wp_idx,:)';
    end

    %% ================= 横向L1 =================
    seg = r_wp_next - r_wp_prev;
    t_hat = seg / max(norm(seg),1e-6);
    t_h = t_hat - dot(t_hat,u_up)*u_up;
    t_h = t_h / max(norm(t_h),1e-6);

    r_from_start = r - r_wp_prev;
    r_xt = r_from_start - dot(r_from_start,t_h)*t_h;
    xt = norm(r_xt);

    if xt > 1e-6
        n_xt = r_xt/xt;
    else
        n_xt = cross(u_up,t_h); n_xt = n_xt/max(norm(n_xt),1e-6);
    end

    % 末段增强：缩L1, 增a_lat_max
    if R_to_T < R_term2
        L1_dist = max(8e3, 8*Vh);
        a_lat_max = a_lat_max_term;
    elseif R_to_T < R_term1
        L1_dist = max(20e3, 12*Vh);
        a_lat_max = 0.8*a_lat_max_term;
    else
        L1_dist = max(L1_base, L1_gainV*Vh);
        a_lat_max = a_lat_max_cruise;
    end

    r_L1 = r + L1_dist*t_h - min(xt,0.5*L1_dist)*n_xt;
    u_L1 = r_L1 - r;
    u_L1 = u_L1 - dot(u_L1,u_up)*u_up;
    u_L1 = u_L1 / max(norm(u_L1),1e-6);

    sin_eta = dot(cross(u_vh,u_L1),u_up);
    cos_eta = dot(u_vh,u_L1);
    eta = atan2(sin_eta,cos_eta);

    a_lat_cmd = 2*Vh^2/max(L1_dist,1)*sin(eta);
    a_lat_cmd = min(max(a_lat_cmd,-a_lat_max),a_lat_max);

    right_h = cross(u_up,u_vh); right_h = right_h/max(norm(right_h),1e-6);
    a_lat_vec = a_lat_cmd * right_h;

    %% ================= 纵向高度PID =================
    [~,~,h_wp_next] = ecef2lla_cgcs2000(r_wp_next');
    h_cmd = h_wp_next;

    v_up = dot(v,u_up);
    h_err = h_cmd - h;

    int_h = int_h + h_err*dt_use;
    int_h = min(max(int_h,-int_h_lim),int_h_lim);

    a_h_cmd = Kph*h_err + Kih*int_h - Kdh*v_up;
    if R_to_T < R_term2
        a_h_lim_now = 1.2*a_h_max;
    else
        a_h_lim_now = a_h_max;
    end
    a_h_cmd = min(max(a_h_cmd,-a_h_lim_now),a_h_lim_now);

    a_vert_vec = a_h_cmd * u_up;
    a_cmd_ecef = a_lat_vec + a_vert_vec;
    a_cmd_norm = norm(a_cmd_ecef);

    %% ================= a_cmd -> alpha/phi =================
    CL_max = 2.8;  % 末段可略提
    aL_cap = (CL_max*qbar*S_ref)/m;
    if a_cmd_norm > aL_cap && a_cmd_norm > 1e-6
        a_cmd_ecef = a_cmd_ecef * (aL_cap/a_cmd_norm);
        a_cmd_norm = aL_cap;
    end

    a_vert_req = dot(a_cmd_ecef,u_up);
    a_h_req_vec = a_cmd_ecef - a_vert_req*u_up;
    a_lat_req = dot(a_h_req_vec,right_h);

    phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
    phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

    CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);

    % 线性反推alpha
    adeg_est = 2.0;
    CL_alpha_rad = (0.07235 - 0.003368*Ma) * (180/pi);
    CL_alpha_rad = max(CL_alpha_rad,0.1);
    CL0_est = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.001185*adeg_est^2;

    alpha_cmd = (CL_req - CL0_est)/CL_alpha_rad;
    alpha_cmd = min(max(alpha_cmd,alpha_min),alpha_max);

    gamma_now = atan2(v_up, Vh);
    theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
    theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);

    chi_v = atan2(dot(u_vh,eE), dot(u_vh,eN));
    psi_cmd = chi_v;

    phi = phi_cmd; theta = theta_cmd; psi = psi_cmd; 

    %% ================= 速度控制（CT前馈+PI） =================
    if R_to_T > 900e3
        M_cmd = M_cmd_far;
    elseif R_to_T > 250e3
        M_cmd = M_cmd_mid;
    else
        M_cmd = M_cmd_near;
    end
    V_cmd = max(M_cmd*a_snd, V_floor);

    [CL,CD,CY,~,~,~] = aero_coeffs(Ma, alpha_cmd, 0, de1, de2, dr, dT);
    CY = 0;

    L = CL*qbar*S_ref;
    D = CD*qbar*S_ref;

    CT_ff = D / max(qbar*S_ref,1);

    v_err = V_cmd - V;
    int_v = int_v + v_err*dt_use;
    int_v = min(max(int_v,-int_v_lim),int_v_lim);

    CT_cmd = CT_ff + Kpv*v_err + Kiv*int_v;
    CT_cmd = min(max(CT_cmd,CT_min),CT_max);

    if (CT_cmd>=CT_max-1e-6 && v_err>0) || (CT_cmd<=CT_min+1e-6 && v_err<0)
        int_v = int_v - 0.7*v_err*dt_use; % anti-windup
    end

    % CT->dT 反解
    dT_grid = linspace(dT_min,dT_max,41);
    CT_grid = zeros(size(dT_grid));
    for ii=1:numel(dT_grid)
        CT_grid(ii) = thrust_coeffs(Ma, alpha_cmd, 0, dT_grid(ii));
    end
    [~,ix] = min(abs(CT_grid - CT_cmd));
    dT_target = dT_grid(ix);

    dT = dT + (dT_target - dT)*dt_use/tau_dT;
    dT = min(max(dT,dT_min),dT_max);

    CT = thrust_coeffs(Ma, alpha_cmd, 0, dT);
    T_eng = CT*qbar*S_ref;

    %% ================= 力与积分 =================
    if norm(v_h) > 1e-6
        fwd = v_h / norm(v_h);
    else
        fwd = v / max(norm(v),1e-6);
    end

    F_drag_e   = -D * (v/max(norm(v),1e-6));
    F_lift_e   =  L * u_up;
    F_thrust_e =  T_eng * fwd;
    F_ecef = F_drag_e + F_lift_e + F_thrust_e;

    if use_simple_gravity
        g_ecef = -mu/norm(r)^3 * r; 
    end

    if use_rotation_terms
        a_ecef = g_ecef + F_ecef/m - 2*cross(omega_ie,v) - cross(omega_ie,cross(omega_ie,r));
    else
        a_ecef = g_ecef + F_ecef/m;
    end

    v = v + a_ecef*dt_use;
    r = r + v*dt_use;
    t_now = t_now + dt_use;

    %% ================= 记录 =================
    t_log(k) = t_now;
    Rhist(k,:) = r.';
    Vhist(k,:) = v.';
    Hhist(k) = h;
    Dist_hist(k) = R_to_T;
    Vc_hist(k) = Vc_toward;
    Vmag_hist(k) = V;
    Vcmd_hist(k) = V_cmd;
    Mcmd_hist(k) = M_cmd;
    dT_hist(k) = dT;
    CL_hist(k) = CL;
    CD_hist(k) = CD;
    CTff_hist(k) = CT_ff;
    CTcmd_hist(k) = CT_cmd;
    L1_hist(k) = L1_dist;
    alat_hist(k) = a_lat_cmd;
    wp_idx_hist(k) = wp_idx;

    if any(~isfinite([r;v;Ma;qbar;CL;CD;CT;dT]))
        stop_reason = "NaN/Inf";
        fprintf('NaN/Inf at t=%.2f s\n', t_now);
        break;
    end

    k = k + 1;
end

%% ===================== 截断有效数据 =====================
valid = isfinite(t_log) & isfinite(Dist_hist);
tt = t_log(valid);

Rk = Rhist(valid,:);
Vk = Vhist(valid,:);
Hk = Hhist(valid);
Dk = Dist_hist(valid);
Vck = Vc_hist(valid);
Vmk = Vmag_hist(valid);
Vcmdk = Vcmd_hist(valid);
Mcmdk = Mcmd_hist(valid);
dTk = dT_hist(valid);
CLk = CL_hist(valid);
CDk = CD_hist(valid);
CTffk = CTff_hist(valid);
CTcmdk = CTcmd_hist(valid);
L1k = L1_hist(valid);
alatk = alat_hist(valid);

fprintf('Simulation stop reason: %s, t=%.2f s, min range=%.1f m\n', stop_reason, tt(end), min(Dk));

%% ===================== 绘图 =====================
figure;
plot(tt, Dk/1000, 'm', 'LineWidth',1.7); grid on;
xlabel('Time (s)'); ylabel('Distance to Target (km)');
title('Distance to Target (Terminal Version)');

figure;
plot(tt, Vck, 'b', 'LineWidth',1.4); grid on; yline(0,'r--');
xlabel('Time (s)'); ylabel('V_c toward target (m/s)');
title('Closing Speed');

figure;
plot(tt, Vmk, 'k', 'LineWidth',1.4); hold on; grid on;
plot(tt, Vcmdk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (m/s)');
title('Speed Tracking');
legend('|V|','V_{cmd}','Location','best');

figure;
plot(tt, Hk/1000, 'LineWidth',1.4); grid on;
xlabel('Time (s)'); ylabel('Altitude (km)');
title('Altitude');

figure;
plot(tt, dTk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('dT');
title('Throttle');

figure;
plot(tt, CTffk, 'b--', 'LineWidth',1.2); hold on; grid on;
plot(tt, CTcmdk, 'r-', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('C_T');
title('CT Feedforward / CT Command');
legend('CT_{ff}','CT_{cmd}','Location','best');

figure;
plot(tt, CLk, 'b', 'LineWidth',1.2); hold on; grid on;
plot(tt, CDk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Coefficient');
title('C_L / C_D');
legend('C_L','C_D','Location','best');

figure;
plot(tt, L1k/1000, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('L1 distance (km)');
title('L1 Distance Scheduling');

figure;
plot(tt, alatk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('a_{lat,cmd} (m/s^2)');
title('Lateral Acceleration Command');

figure;
plot3(Rk(:,1)/1e3, Rk(:,2)/1e3, Rk(:,3)/1e3, 'b','LineWidth',1.3); hold on; grid on; axis equal;
[xe,ye,ze] = sphere(60);
surf(Re*xe/1e3, Re*ye/1e3, Re*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none');
plot3(rT(1)/1e3, rT(2)/1e3, rT(3)/1e3, 'ro','MarkerFaceColor','r');
xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
title('Trajectory in ECEF (Terminal Version)');
legend('Vehicle','Earth','Target','Location','best');