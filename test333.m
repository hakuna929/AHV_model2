% 依赖函数：
% lla2ecef_cgcs2000, ecef2lla_cgcs2000, atmos_simple, aero_coeffs, gravity_cgcs2000_ecef
% 其中 lla2ecef_cgcs2000 / ecef2lla_cgcs2000 的输入输出角度单位均为"度"
% 本脚本内 phi/theta/psi 使用"弧度" alpha使用角度

clear; clc;

%% ================= 仿真时间 =================
dt_base = 0.01;
T_end   = 2000;
N_max   = ceil(T_end/dt_base) * 5;  % 预留给末段小步长

%% ================= 常量 =================
mu  = 3.986004418e14;
Re  = 6378137.0;
we  = 7.2921150e-5;
omega_ie = [0;0;we];
g0 = 9.80665;

%% ================= 飞行器参数 =================
m = 671.33;
S_ref = 0.2986;

%% ================= 起点/终点设置 =================
h0 = 30e3; lat0 = 19.2; lon0 = 110.5;   % lat/lon 单位：度
r = lla2ecef_cgcs2000(lat0, lon0, h0); r = r(:);

latT = 13.35; lonT = 144.55; hT = 30e3; % 度
rT = lla2ecef_cgcs2000(latT, lonT, hT); rT = rT(:);

[a0,~] = atmos_simple(h0);
V0 = 6.5 * a0;

% 航点
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
wp_switch_dist = 100e3;

%% ================= 初始速度方向设置 =================
u_up0 = r/norm(r);
seg0 = r_wp_next - r_wp_prev;
seg0_h = seg0 - dot(seg0,u_up0)*u_up0;
if norm(seg0_h) < 1e-6
    seg0_h = seg0;
end
u_vh0 = seg0_h / norm(seg0_h);

gamma0 = 0*pi/180;
v = V0 * seg0_h / norm(seg0_h);

%% ================= 控制参数设置 =================
% L1 横向
L1_base  = 100e3;
L1_gainV = 25;
a_lat_max_cruise = 60;
a_lat_max_term   = 140;

% ======= 定高：分段PID（按距离R_to_T） =======
R_pid_far = 900e3;
R_pid_mid = 250e3;

PID_far.Kp = 0.010;   PID_far.Ki = 0.00004;   PID_far.Kd = 0.014;
PID_far.int_lim = 8e4;
PID_far.a_max   = 18;

PID_mid.Kp = 0.007;   PID_mid.Ki = 0.000025;  PID_mid.Kd = 0.018;
PID_mid.int_lim = 6e4;
PID_mid.a_max   = 14;

PID_near.Kp = 0.004;  PID_near.Ki = 0.000012; PID_near.Kd = 0.022;
PID_near.int_lim = 4e4;
PID_near.a_max   = 10;

int_h = 0;

phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2;
alpha_max =  10;

alpha_unload_gain = 0;   % 如需饱和卸载，可改为 0.5~2.0

% 速度控制
Kpv = 0.0018;
Kiv = 0.00025;
int_v = 0;
int_v_lim = 8e3;

CT_min = 0.025;
CT_max = 0.3;
dT_min = 0.01;
dT_max = 1.00;
tau_dT = 0.3;
dT = 0.5;

%% ================= 初始配平 =================
use_trim_init = true;
if use_trim_init
    [a0, rho0] = atmos_simple(h0);
    Ma0   = V0 / max(a0,1e-3);
    qbar0 = 0.5 * rho0 * V0^2;

    R0 = 6371000;   % 地球半径 m
    g = g0 * (R0/(R0+h0))^2;
    CL_req0 = m*g / (qbar0*S_ref);

    alpha_trim_min = -2;
    alpha_trim_max = 10;
    alpha0 = fminbnd(@(a) (aero_coeffs(Ma0, a, dT) - CL_req0)^2, alpha_trim_min, alpha_trim_max);

    [~, CD0, ~] = aero_coeffs(Ma0, alpha0, dT);
    D0 = CD0 * qbar0 * S_ref;
    CT_req0 = D0 / max(qbar0*S_ref,1);

    dT_grid = linspace(dT_min, dT_max, 41);
    CT_grid = zeros(size(dT_grid));
    for ii=1:numel(dT_grid)
        [~,~,CT_grid(ii)] = aero_coeffs(Ma0, alpha0, dT_grid(ii));
    end
    [~,ix] = min(abs(CT_grid - CT_req0));
    dT = dT_grid(ix);

    eE0 = [-sind(lon0); cosd(lon0); 0];
    eN0 = [-sind(lat0)*cosd(lon0); -sind(lat0)*sind(lon0); cosd(lat0)];
    chi0 = atan2(dot(u_vh0,eE0), dot(u_vh0,eN0));

    phi0   = 0;
    theta0 = gamma0 + deg2rad(alpha0);
    psi0   = chi0;

    fprintf('Trim init: alpha0=%.3f deg, dT0=%.3f, theta0=%.3f deg\n', ...
        alpha0, dT, theta0*180/pi);
end

%% ================= 马赫指令 =================
M_cmd_far  = 6.5;
M_cmd_mid  = 5.5;
M_cmd_near = 5.0;
V_floor    = 1100;

%% ================= 终端捕获参数 =================
R_hit            = 1000;
R_term1          = 120e3;
R_term2          = 40e3;
R_pass_check     = 20e3;
R_abort_diverge  = 50e3;

%% ================= 记录数组 =================
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
CT_hist = nan(N_max,1);
CTff_hist = nan(N_max,1);
CTcmd_hist = nan(N_max,1);
L1_hist = nan(N_max,1);
alat_hist = nan(N_max,1);
wp_idx_hist = nan(N_max,1);
alpha_hist = nan(N_max,1);
phi_hist   = nan(N_max,1);
theta_hist = nan(N_max,1);
T_hist     = nan(N_max,1);
D_hist     = nan(N_max,1);
L_hist     = nan(N_max,1);

stop_reason = "completed";
k = 1;
t_now = 0;

prev_R = inf;
min_R = inf;

% 用于解耦CL_req时的滚转补偿：这里用上一时刻phi
phi = 0; theta = 0; psi = 0;

%% ================= 主循环 =================
while k <= N_max && t_now <= T_end

    [lat_deg, lon_deg, h] = ecef2lla_cgcs2000(r');

    if h <= 20e3
        stop_reason = "terrain impact";
        fprintf('Vehicle impacted terrain at t=%.1f s\n', t_now);
        break;
    end

    r_rel_T = rT - r;
    R_to_T = norm(r_rel_T);
    u_los_T = r_rel_T / max(R_to_T,1);
    Vc_toward = dot(v, u_los_T);

    min_R = min(min_R, R_to_T);

    if R_to_T <= R_hit
        stop_reason = "target reached";
        fprintf('Target reached at t=%.1f s, miss=%.1f m\n', t_now, R_to_T);
        break;
    end

    if (R_to_T < R_pass_check) && (Vc_toward <= 0)
        stop_reason = "passed closest approach";
        fprintf('Passed closest approach at t=%.1f s, min range=%.1f m\n', t_now, min_R);
        break;
    end

    if (prev_R < R_abort_diverge) && (R_to_T > prev_R + 30)
        stop_reason = "diverging after near pass";
        fprintf('Diverging after near pass at t=%.1f s, range=%.1f m\n', t_now, R_to_T);
        break;
    end
    prev_R = R_to_T;

    if R_to_T < R_term2
        dt_use = 0.001;
    elseif R_to_T < R_term1
        dt_use = 0.01;
    else
        dt_use = dt_base;
    end

    %% 大气
    V = max(norm(v),1e-3);
    [a_snd, rho] = atmos_simple(max(h,0));
    Ma = max(min(V/max(a_snd,1e-3),8.0),0.0);
    qbar = 0.5*rho*V^2;

    %% 局部基
    u_up = r / norm(r);
    eE = [-sind(lon_deg); cosd(lon_deg); 0];
    eN = [-sind(lat_deg)*cosd(lon_deg); -sind(lat_deg)*sind(lon_deg); cosd(lat_deg)];

    v_h = v - dot(v,u_up)*u_up;
    Vh = max(norm(v_h),1e-6);
    u_vh = v_h / Vh;

    %% 航点切换
    if norm(r_wp_next - r) < wp_switch_dist && wp_idx < nWP
        wp_idx = wp_idx + 1;
        r_wp_prev = r_wp_next;
        r_wp_next = wps_ecef(wp_idx,:)';
    end

    %% ================= 横向 L1 =================
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
        n_xt = cross(u_up,t_h);
        n_xt = n_xt/max(norm(n_xt),1e-6);
    end

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

    right_h = cross(u_up,u_vh);
    right_h = right_h/max(norm(right_h),1e-6);
    a_lat_vec = a_lat_cmd * right_h;

    %% ================= 高度环（分段PID + 实时重力前馈） =================
    % 分段选择
    if R_to_T > R_pid_far
        PID = PID_far;
    elseif R_to_T > R_pid_mid
        PID = PID_mid;
    else
        PID = PID_near;
    end
    Kph = PID.Kp; Kih = PID.Ki; Kdh = PID.Kd;
    int_h_lim = PID.int_lim;
    a_h_max   = PID.a_max;

    h_cmd = 30e3;

    % 实时重力（ECEF）& 沿天顶方向的向下重力大小（正数）
    [lat_g, lon_g, h_g] = ecef2lla_cgcs2000(r');
    g_ecef = gravity_cgcs2000_ecef(lat_g, lon_g, h_g);
    g_up = -dot(g_ecef, u_up);

    v_up = dot(v,u_up);
    h_err = h_cmd - h;

    int_h = int_h + h_err*dt_use;
    int_h = min(max(int_h,-int_h_lim),int_h_lim);

    % a_h_cmd 是"沿u_up方向的加速度指令"（包含重力前馈）
    a_h_cmd = Kph*h_err + Kih*int_h - Kdh*v_up + g_up;

    % 限幅（保留你原来近端放大一点的习惯）
    if R_to_T < R_term2
        a_h_lim_now = 1.2*a_h_max;
    else
        a_h_lim_now = a_h_max;
    end
    a_h_cmd = min(max(a_h_cmd,-a_h_lim_now),a_h_lim_now);

    % 前200s高度保护
    h_floor_200 = 29e3;
    t_protect   = 2;
    if (t_now < t_protect) && (h < h_floor_200)
        a_h_cmd = max(a_h_cmd, 0);
        int_h   = max(int_h, 0);
    end

    % 合成指令加速度（仍保留给后面算phi）
    a_vert_vec = a_h_cmd * u_up;
    a_cmd_ecef = a_lat_vec + a_vert_vec;

    %% ================= 解耦：用"垂向需求"反解 CL_req -> alpha_cmd =================
    % 净向上加速度需求（气动需要提供的部分）
    % a_up_net = a_h_cmd - g_up;
    % 
    % % 可选：限制向下净加速度，避免过度压杆导致耦合
    % a_up_net = max(a_up_net, -2.0);
    % 
    % % 滚转导致垂向升力损失：L_up = L*cos(phi)
    % cphi = cos(abs(phi));      % 用上一时刻phi（本循环phi_cmd还未算）
    % cphi = max(cphi, 0.35);    % 保护，避免除0/过大放大
    % 
    % L_up_req = m * a_up_net;
    % L_req = L_up_req / cphi;
    % 
    % CL_req = L_req / max(qbar*S_ref, 1.0);
    % CL_req = min(max(CL_req, -0.5), 1.2);   % 可按你的气动范围调整

    % a_h_cmd 是你希望沿 u_up 的"总加速度"（如果你在PID里加了g_up前馈，它通常接近g_up）
    a_up_total = a_h_cmd;

    % 防止数值异常：至少保证不小于某个比例的g_up，否则会持续下沉
    a_up_total = max(a_up_total, 0.8*g_up);

    cphi = max(cos(abs(phi)), 0.35);
    L_req  = m * a_up_total / cphi;
    CL_req = L_req / max(qbar*S_ref, 1.0);

    % 用线性近似反解alpha
    % adeg_est = 0.0;
    % CL_alpha_deg = (0.07235 - 0.003368*Ma) * (pi/180);
    % CL_alpha_deg = max(CL_alpha_deg,0.1);
    % CL0_est = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.001185*adeg_est^2;

    % [alpha_cmd, invInfo] = invert_alpha_from_CL(Ma, dT, CL_req, alpha_min, alpha_max);

    alpha_cmd = fminbnd(@(a) (aero_coeffs(Ma, a, dT) - CL_req).^2, alpha_min, alpha_max);

    % 
    % alpha_cmd = (CL_req - CL0_est)/CL_alpha_deg;
    alpha_sat = min(max(alpha_cmd, alpha_min), alpha_max);

    % 饱和卸载（可选）
    if abs(alpha_cmd - alpha_sat) > 1e-9
        over = abs(alpha_cmd - alpha_sat)/max(abs(alpha_max-alpha_min),1e-6);
        unload = min(1.0, alpha_unload_gain * (0.2 + 3.0*over));
        a_h_cmd = (1.0 - unload)*a_h_cmd;
        int_h = int_h * (1.0 - 0.5*unload);
        alpha_cmd = alpha_sat;
    else
        alpha_cmd = alpha_sat;
    end

    %% ================= a_cmd -> phi/theta/psi（保持你原结构） =================
    a_vert_req = dot(a_cmd_ecef,u_up);
    a_h_req_vec = a_cmd_ecef - a_vert_req*u_up;
    a_lat_req = dot(a_h_req_vec,right_h);

    phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
    phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

    gamma_now = atan2(v_up, Vh);
    theta_cmd = gamma_now + deg2rad(alpha_cmd)*cos(phi_cmd);
    theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);

    chi_v = atan2(dot(u_vh,eE), dot(u_vh,eN));
    psi_cmd = chi_v;

    phi = phi_cmd;
    theta = theta_cmd;
    psi = psi_cmd;

    %% ================= 速度控制（CT 前馈 + PI） =================
    if R_to_T > 900e3
        M_cmd = M_cmd_far;
    elseif R_to_T > 250e3
        M_cmd = M_cmd_mid;
    else
        M_cmd = M_cmd_near;
    end
    V_cmd = max(M_cmd*a_snd, V_floor);

    [CL, CD, ~] = aero_coeffs(Ma, alpha_cmd, dT);
    L = CL*qbar*S_ref;
    D = CD*qbar*S_ref;

    CT_ff = D / max(qbar*S_ref,1);

    v_err = V_cmd - V;
    int_v = int_v + v_err*dt_use;
    int_v = min(max(int_v,-int_v_lim),int_v_lim);

    CT_cmd = CT_ff + Kpv*v_err + Kiv*int_v;
    CT_cmd = min(max(CT_cmd,CT_min),CT_max);

    if (CT_cmd>=CT_max-1e-6 && v_err>0) || (CT_cmd<=CT_min+1e-6 && v_err<0)
        int_v = int_v - 0.7*v_err*dt_use;
    end

    dT_grid = linspace(dT_min,dT_max,41);
    CT_grid = zeros(size(dT_grid));
    for ii=1:numel(dT_grid)
        [~,~,CT_grid(ii)] = aero_coeffs(Ma, alpha_cmd, dT_grid(ii));
    end
    [~,ix] = min(abs(CT_grid - CT_cmd));
    dT_target = dT_grid(ix);

    dT = dT + (dT_target - dT)*dt_use/tau_dT;
    dT = min(max(dT,dT_min),dT_max);

    [~,~,CT] = aero_coeffs(Ma, alpha_cmd, dT);
    T_eng = CT * qbar * S_ref;

    %% ================= 力与积分（3D 升力方向含滚转phi） =================
    Vmag = max(norm(v), 1e-6);
    u_v  = v / Vmag;

    right = cross(u_up, u_v);
    nr = norm(right);
    if nr < 1e-8
        right = cross(u_up, eE);
        nr = norm(right);
    end
    right = right / max(nr, 1e-6);

    lift_up = cross(u_v, right);
    lift_up = lift_up / max(norm(lift_up), 1e-6);

    lift_dir = cos(phi) * lift_up + sin(phi) * right;
    lift_dir = lift_dir / max(norm(lift_dir), 1e-6);

    F_drag_e   = -D * u_v;
    F_lift_e   =  L * lift_dir;
    F_thrust_e =  T_eng * u_v;

    F_ecef = F_drag_e + F_lift_e + F_thrust_e;

    % 复用上面高度环算的g_ecef（保持一致）
    a_ecef = g_ecef + F_ecef/m ...
        - 2*cross(omega_ie, v) ...
        - cross(omega_ie, cross(omega_ie, r));

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
    CT_hist(k) = CT;
    CTff_hist(k) = CT_ff;
    CTcmd_hist(k) = CT_cmd;
    L1_hist(k) = L1_dist;
    alat_hist(k) = a_lat_cmd;
    wp_idx_hist(k) = wp_idx;
    alpha_hist(k) = alpha_cmd;
    phi_hist(k)   = phi_cmd;
    theta_hist(k) = theta_cmd;
    T_hist(k)     = T_eng;
    D_hist(k)     = D;
    L_hist(k)     = L;
    % ===== 新增：升力/重力对比日志 =====
    Lmag_hist   = nan(N_max,1);   % |L| (N)
    Wmag_hist   = nan(N_max,1);   % m*g_up (N)
    Lup_hist    = nan(N_max,1);   % 升力在u_up方向分量 (N)
    Wup_hist    = nan(N_max,1);   % 重力在u_up方向分量大小(向下为正) (N)
    CLreq_hist  = nan(N_max,1);   % 你用于反解/控制的CL_req

    % ===== 新增记录：升力/重力对比 =====
    CLreq_hist(k) = CL_req;          % 你当前用的CL_req（无论线性/插值/优化反解都一样记）
    Lmag_hist(k)  = abs(L);          % 升力大小（N）
    Wmag_hist(k)  = m * g_up;        % 重力大小（N）

    % 注意：你这里的"升力力向量"是 F_lift_e = L * lift_dir
    % 所以其在天顶u_up方向的分量：
    Lup_hist(k) = dot(F_lift_e, u_up);

    % 重力在u_up方向的分量（向下为正）：g_ecef·u_up通常为负，所以取负号
    Wup_hist(k) = -m * dot(g_ecef, u_up);

    k = k + 1;
end

%% ================= 截断有效数据 =================
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
Lmagk  = Lmag_hist(valid);
Wmagk  = Wmag_hist(valid);
Lupk   = Lup_hist(valid);
Wupk   = Wup_hist(valid);
CLreqk = CLreq_hist(valid);
CDk = CD_hist(valid);
CTk = CT_hist(valid);
CTffk = CTff_hist(valid);
CTcmdk = CTcmd_hist(valid);
L1k = L1_hist(valid);
alatk = alat_hist(valid);
alphak = alpha_hist(valid);
phik   = phi_hist(valid);
thetak = theta_hist(valid);
Tk     = T_hist(valid);
Dforce = D_hist(valid);
Lforce = L_hist(valid);

fprintf('Simulation stop reason: %s, t=%.2f s, min range=%.1f m\n', stop_reason, tt(end), min(Dk));

%% ================= 绘图 =================
figure;
plot(tt, Vmk, 'k', 'LineWidth',1.4); hold on; grid on;
plot(tt, Vcmdk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (m/s)');
title('速度跟踪');
legend('|V|','V_{cmd}','Location','best');

figure;
plot(tt, Hk/1000, 'LineWidth',1.4); grid on; yline(30,'r--');
xlabel('Time (s)'); ylabel('Altitude (km)');
title('高度保持');

figure;
plot(tt, alphak, 'LineWidth',1.3); grid on;
yline(alpha_min,'r--'); yline(alpha_max,'r--');
xlabel('Time (s)'); ylabel('\alpha (deg)');
title('攻角');

figure;
plot(tt, dTk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('dT');
title('油门开度');

figure;
plot(tt, CTffk, 'b--', 'LineWidth',1.2); hold on; grid on;
plot(tt, CTcmdk, 'r-', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('C_T');
title('CT 前馈 / CT 指令');
legend('CT_{ff}','CT_{cmd}','Location','best');

figure;
plot(tt, CTk, 'k', 'LineWidth',1.2); grid on;
xlabel('Time (s)'); ylabel('C_T (actual)');
title('CT Actual');

figure;
plot(tt, CLk, 'b', 'LineWidth',1.3); hold on; grid on;
plot(tt, CLreqk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('C_L');
title('升力系数：C_L(实际) vs C_{L,req}');
legend('C_L','C_{L,req}','Location','best');

figure;
plot(tt, Lmagk/1e3, 'b', 'LineWidth',1.3); hold on; grid on;
plot(tt, Wmagk/1e3, 'r--', 'LineWidth',1.3);
xlabel('Time (s)'); ylabel('Force (kN)');
title('升力/重力大小对比');
legend('|L|','m g','Location','best');

figure;
plot(tt, Lupk/1e3, 'b', 'LineWidth',1.3); hold on; grid on;
plot(tt, Wupk/1e3, 'r--', 'LineWidth',1.3);
xlabel('Time (s)'); ylabel('Force along u_{up} (kN)');
title('天顶方向力平衡：L_{up} vs W_{up}');
legend('L_{up}','W_{up}','Location','best');

figure;
plot(tt, CLk, 'b', 'LineWidth',1.2); hold on; grid on;
plot(tt, CDk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Coefficient');
title('升阻比');
legend('C_L','C_D','Location','best');

figure;
plot(tt, alatk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('a_{lat,cmd} (m/s^2)');
title('侧向加速度指令');

figure;
plot(tt, Tk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('Thrust (N)');
title('发动机推力');

%% ================= 绘图（发射系 ENU） =================
r0 = lla2ecef_cgcs2000(lat0, lon0, h0).';
r0 = r0(:);

slat = sind(lat0); clat = cosd(lat0);
slon = sind(lon0); clon = cosd(lon0);

C_ecef2enu = [ -slon,        clon,       0;
    -slat*clon,  -slat*slon,   clat;
    clat*clon,   clat*slon,   slat ];

dR = (Rk.' - r0);
enu = C_ecef2enu * dR;

E = enu(1,:)/1e3;
N = enu(2,:)/1e3;
U = enu(3,:)/1e3;

figure;
plot3(E, N, U, 'b', 'LineWidth', 1.5); grid on; axis equal;
xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
title('发射系下轨迹 (ENU)');

%% ================= 轨迹经纬度图 =================
Ntraj = size(Rk,1);
lat_traj = nan(Ntraj,1);
lon_traj = nan(Ntraj,1);
h_traj   = nan(Ntraj,1);

for ii = 1:Ntraj
    [lat_traj(ii), lon_traj(ii), h_traj(ii)] = ecef2lla_cgcs2000(Rk(ii,:));
end

figure;
plot3(lon_traj, lat_traj, h_traj/1000, 'b', 'LineWidth', 1.5); grid on; hold on;
plot3(lon0, lat0, h0/1000, 'go', 'MarkerFaceColor','g');
plot3(lonT, latT, hT/1000, 'ro', 'MarkerFaceColor','r');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); zlabel('Altitude (km)');
title('经纬高轨迹 LLA (Lon-Lat-Alt)');
legend('Trajectory','Launch','Target','Location','best');
view(3);
