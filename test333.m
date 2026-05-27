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

% 定高
Kph = 0.010;
Kih = 0.00005;
Kdh = 0.012;
int_h = 0;
int_h_lim = 8e4;
a_h_max = 18;

use_gravity_ff = true;
a_h_bias = 1.2*g0;

phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2;
alpha_max =  10;

alpha_unload_gain = 0.1;

% 速度控制
Kpv = 0.0018;
Kiv = 0.00025;
int_v = 0;
int_v_lim = 8e3;

CT_min = 0.025;
CT_max = 0.09;
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
    % CL_req0 = m*g0 / max(qbar0*S_ref,1);

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
M_cmd_mid  = 6.5;
M_cmd_near = 6.5;
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
m_hist     = nan(N_max,1);
mf_hist    = nan(N_max,1);
mdot_hist  = nan(N_max,1);

stop_reason = "completed";
k = 1;
t_now = 0;

prev_R = inf;
min_R = inf;

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

    %% ================= 定高（纵向） =================
    h_cmd = 30e3;

    v_up = dot(v,u_up);
    h_err = h_cmd - h;

    int_h = int_h + h_err*dt_use;
    int_h = min(max(int_h,-int_h_lim),int_h_lim);

    if use_gravity_ff
        a_h_cmd = Kph*h_err + Kih*int_h - Kdh*v_up + a_h_bias;
    else
        a_h_cmd = Kph*h_err + Kih*int_h - Kdh*v_up;
    end

    if R_to_T < R_term2
        a_h_lim_now = 1.2*a_h_max;
    else
        a_h_lim_now = a_h_max;
    end
    a_h_cmd = min(max(a_h_cmd,-a_h_lim_now),a_h_lim_now);

    %% 前200s高度保护
    h_floor_200 = 29e3;
    t_protect   = 2;
    if (t_now < t_protect) && (h < h_floor_200)
        a_h_cmd = max(a_h_cmd, 0);
        int_h   = max(int_h, 0);
    end

    a_vert_vec = a_h_cmd * u_up;
    a_cmd_ecef = a_lat_vec + a_vert_vec;
    a_cmd_norm = norm(a_cmd_ecef);

    %% ================= a_cmd -> alpha/phi =================
    CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);

    adeg_est = 0.0;
    CL_alpha_deg = (0.07235 - 0.003368*Ma) * (pi/180);
    CL_alpha_deg = max(CL_alpha_deg,0.1);
    CL0_est = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.001185*adeg_est^2;

    alpha_cmd = (CL_req - CL0_est)/CL_alpha_deg;
    alpha_sat = min(max(alpha_cmd, alpha_min), alpha_max);

    if abs(alpha_cmd - alpha_sat) > 1e-9
        over = abs(alpha_cmd - alpha_sat)/max(abs(alpha_max-alpha_min),1e-6);
        unload = min(1.0, alpha_unload_gain * (0.2 + 3.0*over));
        a_h_cmd = (1.0 - unload)*a_h_cmd;
        int_h = int_h * (1.0 - 0.5*unload);

        a_vert_vec = a_h_cmd * u_up;
        a_cmd_ecef = a_lat_vec + a_vert_vec;
        a_cmd_norm = norm(a_cmd_ecef);

        CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);
        alpha_cmd = (CL_req - CL0_est)/CL_alpha_deg;
        alpha_cmd = min(max(alpha_cmd, alpha_min), alpha_max);
    else
        alpha_cmd = alpha_sat;
    end

    a_vert_req = dot(a_cmd_ecef,u_up);
    a_h_req_vec = a_cmd_ecef - a_vert_req*u_up;
    a_lat_req = dot(a_h_req_vec,right_h);

    phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
    phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

    gamma_now = atan2(v_up, Vh);
    theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
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

    %% ================= 质量更新（正式调用独立函数） =================
    % [m, m_fuel, mdot_f] = update_mass_model( ...
    %     m_struct, m_fuel, Ma, alpha_cmd, 0.0, h/1000, dT, dt_use, ...
    %     mass_model_alpha_in_deg, mdot_min, mdot_max);
    % 
    % if m <= m_struct
    %     m = m_struct;
    %     m_fuel = 0;
    %     mdot_f = 0;
    % end

    %% ================= 力与积分 =================
    if norm(v_h) > 1e-6
        fwd = v_h / norm(v_h);
    else
        fwd = v / max(norm(v),1e-6);
    end

    F_drag_e   = -D * (v/max(norm(v),1e-6));
    F_lift_e   =  L * u_up;
    F_thrust_e  =  T_eng * fwd;
    F_ecef = F_drag_e + F_lift_e + F_thrust_e;

    [lat_g, lon_g, h_g] = ecef2lla_cgcs2000(r');
    g_ecef = gravity_cgcs2000_ecef(lat_g, lon_g, h_g);

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
    % m_hist(k)     = m;
    % mf_hist(k)    = m_fuel;
    % mdot_hist(k)  = mdot_f;

    % if any(~isfinite([r;v;Ma;qbar;CL;CD;CT;dT;m;m_fuel;mdot_f]))
    %     stop_reason = "NaN/Inf";
    %     fprintf('NaN/Inf at t=%.2f s\n', t_now);
    %     break;
    % end

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
% mk     = m_hist(valid);
% mfk    = mf_hist(valid);
% mdotk  = mdot_hist(valid);

fprintf('Simulation stop reason: %s, t=%.2f s, min range=%.1f m\n', stop_reason, tt(end), min(Dk));

%% ================= 绘图 =================
% figure;
% plot(tt, Dk/1000, 'm', 'LineWidth',1.7); grid on;
% xlabel('Time (s)'); ylabel('Distance to Target (km)');
% title('Distance to Target');

% figure;
% plot(tt, Vck, 'b', 'LineWidth',1.4); grid on; yline(0,'r--');
% xlabel('Time (s)'); ylabel('V_c toward target (m/s)');
% title('Closing Speed');

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
plot(tt, CLk, 'b', 'LineWidth',1.2); hold on; grid on;
plot(tt, CDk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Coefficient');
title('升阻比');
legend('C_L','C_D','Location','best');

figure;
plot(tt, L1k/1000, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('L1 distance (km)');
title('L1距离');

figure;
plot(tt, alatk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('a_{lat,cmd} (m/s^2)');
title('侧向加速度指令');

figure;
plot(tt, Tk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('Thrust (N)');
title('发动机推力');

% figure;
% plot(tt, (Tk - Dforce), 'LineWidth',1.3); grid on;
% xlabel('Time (s)'); ylabel('(T - D) (N)');
% title('Excess Thrust');

% figure;
% plot(tt, mk, 'LineWidth',1.3); grid on;
% xlabel('Time (s)'); ylabel('Mass (kg)');
% title('Total Mass');

% figure;
% plot(tt, mfk, 'LineWidth',1.3); grid on;
% xlabel('Time (s)'); ylabel('Fuel Mass (kg)');
% title('Fuel Mass');

% figure;
% plot(tt, mdotk, 'LineWidth',1.3); grid on;
% xlabel('Time (s)'); ylabel('Fuel Mass Flow (kg/s)');
% title('燃料消耗');

% figure;
% plot3(Rk(:,1)/1e3, Rk(:,2)/1e3, Rk(:,3)/1e3, 'b','LineWidth',1.3); hold on; grid on; axis equal;
% [xe,ye,ze] = sphere(60);
% surf(Re*xe/1e3, Re*ye/1e3, Re*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none');
% plot3(rT(1)/1e3, rT(2)/1e3, rT(3)/1e3, 'ro','MarkerFaceColor','r');
% xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
% title('Trajectory in ECEF');
% legend('Vehicle','Earth','Target','Location','best');

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

% figure;
% plot(lon_traj, lat_traj, 'b', 'LineWidth', 1.5); grid on; hold on;
% plot(lon0, lat0, 'go', 'MarkerFaceColor','g');
% plot(lonT, latT, 'ro', 'MarkerFaceColor','r');
% xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
% title('Trajectory Ground Track (Lat/Lon)');
% legend('Trajectory','Launch','Target','Location','best');

figure;
plot3(lon_traj, lat_traj, h_traj/1000, 'b', 'LineWidth', 1.5); grid on; hold on;
plot3(lon0, lat0, h0/1000, 'go', 'MarkerFaceColor','g');
plot3(lonT, latT, hT/1000, 'ro', 'MarkerFaceColor','r');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); zlabel('Altitude (km)');
title('经纬高轨迹 LLA (Lon-Lat-Alt)');
legend('Trajectory','Launch','Target','Location','best');
view(3);

% figure;
% subplot(2,1,1);
% plot(lon_traj, h_traj/1000, 'LineWidth', 1.3); grid on;
% xlabel('Longitude (deg)'); ylabel('Altitude (km)'); title('Altitude vs Longitude');
% 
% subplot(2,1,2);
% plot(lat_traj, h_traj/1000, 'LineWidth', 1.3); grid on;
% xlabel('Latitude (deg)'); ylabel('Altitude (km)'); title('Altitude vs Latitude');

%% ================= local functions =================
function CL = getCL(Ma, alpha, dT)
    [CL,~,~] = aero_coeffs(Ma, alpha, dT);
end