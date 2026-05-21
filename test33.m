% main_3dof_MPC_CL_hold.m
% 修正版：
% - 横向：轻量MPC（离散候选集）
% - 纵向：改为"CL需求保持"，由高度误差直接生成CL_req，再反推alpha
% - 保留速度控制、推力反解、攻角/姿态限幅、终端判据
%
% 目标：
% 1) 避免一直掉高
% 2) 避免 theta_cmd 变成向量
% 3) 避免 MPC 和高度环互相打架
%
% 依赖函数：
% lla2ecef_cgcs2000, ecef2lla_cgcs2000, atmos_simple, aero_coeffs

clear; clc;

%% ========== 仿真时间 ==========
dt_base = 0.01;
T_end   = 3000;
N_max   = ceil(T_end/dt_base) * 5;

%% ========== 常量 ==========
mu  = 3.986004418e14;
Re  = 6378137.0;
we  = 7.2921150e-5;
omega_ie = [0;0;we];
g0 = 9.80665;

%% ========== 飞行器参数 ==========
m     = 671.33;
S_ref = 0.2986;

%% ========== 起点/终点 ==========
h0 = 25e3; lat0 = 19.2; lon0 = 110.5;
r = lla2ecef_cgcs2000(lat0, lon0, h0); r = r(:);

latT = 13.35; lonT = 144.55; hT = 25e3;
rT = lla2ecef_cgcs2000(latT, lonT, hT); rT = rT(:);

[a0,~] = atmos_simple(h0);
V0 = 6.5 * a0;

%% 初始速度方向
u_up0 = r/norm(r);
seg0 = rT - r;
seg0_h = seg0 - dot(seg0,u_up0)*u_up0;
if norm(seg0_h) < 1e-6, seg0_h = seg0; end
u_vh0 = seg0_h / norm(seg0_h);

gamma0 = 0*pi/180;
v = V0 * seg0_h / norm(seg0_h);

%% ========== 横向轻量 MPC ==========
a_lat_candidates = [-60, -40, -20, 0, 20, 40, 60];
a_lat_max_cruise = 60;
a_lat_max_term   = 140;

%% ========== 纵向：CL保持 ==========
h_cmd_cruise = 25e3;

% CL需求控制增益（比原来的 a_h 命令更稳）
Kp_CL = 2.0e-5;
Ki_CL = 1.0e-8;
Kd_CL = 5.0e-6;

int_h = 0;
int_h_lim = 8e4;

% 给CL需求加一个基础值，接近平飞配平
CL_floor = 0.15;
CL_ceil  = 2.0;

% alpha/phi/theta限幅
phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2*pi/180;
alpha_max =  10*pi/180;
alpha_unload_gain = 0.1;

%% ========== 速度控制 ==========
M_cmd_far  = 6;
M_cmd_mid  = 5.5;
M_cmd_near = 5;
V_floor    = 1100;

Kpv = 0.0018;
Kiv = 0.00025;
int_v = 0;
int_v_lim = 8e3;

CT_min = 0.025; CT_max = 0.09;
dT_min = 0.001; dT_max = 1.00;
tau_dT = 0.3; dT = 0.5;

%% ========== 终端判据 ==========
R_hit            = 1000;
R_term1          = 120e3;
R_term2          = 40e3;
R_pass_check     = 20e3;
R_abort_diverge  = 50e3;

%% ========== 记录数组 ==========
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
aLat_hist = nan(N_max,1);
CLreq_hist = nan(N_max,1);
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

%% ========== 初始配平 ==========
use_trim_init = true;
if use_trim_init
    [a0, rho0] = atmos_simple(h0);
    Ma0   = V0 / max(a0,1e-3);
    qbar0 = 0.5 * rho0 * V0^2;

    CL_req0 = m*g0 / max(qbar0*S_ref,1);

    alpha_trim_min = -5*pi/180;
    alpha_trim_max = 10*pi/180;
    alpha0 = fminbnd(@(a) (getCL(Ma0, a, dT) - CL_req0).^2, ...
                     alpha_trim_min, alpha_trim_max);

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
    theta0 = gamma0 + alpha0;
    psi0   = chi0;

    fprintf('Trim init: alpha0=%.3f deg, dT0=%.3f, theta0=%.3f deg\n', ...
            alpha0*180/pi, dT, theta0*180/pi);
end

%% ========== 可行性检查 ==========
[a_chk, rho_chk] = atmos_simple(h0);
Ma_chk = V0 / max(a_chk,1e-3);
qbar_chk = 0.5 * rho_chk * V0^2;

CL_req_chk = m * g0 / max(qbar_chk * S_ref, 1);

% 取最大攻角下的可用升力
[CL_max_chk, ~, ~] = aero_coeffs(Ma_chk, alpha_max, dT);

% 取最小攻角下的可用升力（用于看范围）
[CL_min_chk, ~, ~] = aero_coeffs(Ma_chk, alpha_min, dT);

fprintf('\n========== Feasibility Check ==========\n');
fprintf('Initial Mach       : %.3f\n', Ma_chk);
fprintf('Required CL for 30km: %.3f\n', CL_req_chk);
fprintf('Available CL @ alpha_min(%.1f deg): %.3f\n', alpha_min*180/pi, CL_min_chk);
fprintf('Available CL @ alpha_max(%.1f deg): %.3f\n', alpha_max*180/pi, CL_max_chk);

if CL_max_chk < CL_req_chk
    fprintf('>>> INFEASIBLE: alpha_max is not enough to hold 30 km altitude.\n');
    fprintf('>>> Suggestion: increase alpha_max or reduce h_cmd_cruise.\n');

    % 估算一个可行高度（粗略）
    h_test = h0;
    for hh = 20e3:500:40e3
        [a_hh, rho_hh] = atmos_simple(hh);
        q_hh = 0.5 * rho_hh * V0^2;
        CL_req_hh = m * g0 / max(q_hh * S_ref,1);
        [CL_max_hh,~,~] = aero_coeffs(Ma_chk, alpha_max, dT);
        if CL_max_hh >= CL_req_hh
            h_test = hh;
            break;
        end
    end
    fprintf('Approx feasible altitude with alpha_max=%.1f deg: %.1f km\n', ...
            alpha_max*180/pi, h_test/1000);
else
    fprintf('>>> FEASIBLE: max alpha can hold the initial altitude.\n');
end
fprintf('======================================\n\n');

%% ===================== 主循环 =====================
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

    V = max(norm(v),1e-3);
    [a_snd, rho] = atmos_simple(max(h,0));
    Ma = max(min(V/max(a_snd,1e-3),8.0),0.0);
    qbar = 0.5*rho*V^2;

    u_up = r / norm(r);
    eE = [-sind(lon_deg); cosd(lon_deg); 0];
    eN = [-sind(lat_deg)*cosd(lon_deg); -sind(lat_deg)*sind(lon_deg); cosd(lat_deg)];

    v_up = dot(v,u_up);
    v_h_vec = v - v_up*u_up;
    Vh = max(norm(v_h_vec),1e-6);
    u_vh = v_h_vec / Vh;

    %% ========== 横向轻量 MPC ==========
    if R_to_T < R_term2
        a_lat_max = a_lat_max_term;
    elseif R_to_T < R_term1
        a_lat_max = 0.8*a_lat_max_term;
    else
        a_lat_max = a_lat_max_cruise;
    end

    a_lat_cand = a_lat_candidates(abs(a_lat_candidates) <= a_lat_max);
    if isempty(a_lat_cand), a_lat_cand = 0; end

    right_h = cross(u_up,u_vh);
    right_h = right_h / max(norm(right_h),1e-6);

    bestJ = inf;
    best_a_lat = 0;
    best_alpha = 0;
    best_phi = 0;
    best_theta = 0;
    best_CL = 0;
    best_CD = 0;
    best_CT = 0;
    best_CTff = 0;
    best_CLreq = 0;
    best_dT = dT;

    %% ========== 纵向：CL需求 ==========
    h_err = h_cmd_cruise - h;
    int_h = int_h + h_err*dt_use;
    int_h = min(max(int_h,-int_h_lim),int_h_lim);

    % 基于高度误差直接生成CL需求，避免先算a_h_cmd再反推
    CL_trim_req = m*g0 / max(qbar*S_ref,1);
    CL_req_base = CL_trim_req + Kp_CL*h_err + Ki_CL*int_h - Kd_CL*v_up;

    % 前200s高度保护
    if t_now < 2 && h < 29e3
        CL_req_base = max(CL_req_base, CL_trim_req);
        int_h = max(int_h, 0);
    end

    CL_req_base = min(max(CL_req_base, CL_floor), CL_ceil);

    gamma_now = atan2(v_up, Vh);
    if ~isscalar(gamma_now), gamma_now = gamma_now(1); end

    for ia = 1:numel(a_lat_cand)

        a_lat_cmd = a_lat_cand(ia);

        a_cmd_norm_guess = sqrt(a_lat_cmd^2 + (CL_req_base*qbar*S_ref/m)^2);
        CL_req = (m*a_cmd_norm_guess)/max(qbar*S_ref,1);

        % 把纵向和横向合在一起估算所需CL
        a_vert_guess = CL_req_base * qbar*S_ref/m;
        a_cmd_ecef_guess = a_lat_cmd*right_h + a_vert_guess*u_up;
        a_cmd_norm = norm(a_cmd_ecef_guess);

        CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);

        % 线性反推alpha
        adeg_est = 0.0;
        CL_alpha_rad = (0.07235 - 0.003368*Ma) * (180/pi);
        CL_alpha_rad = max(CL_alpha_rad,0.1);
        CL0_est = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.001185*adeg_est^2;

        alpha_raw = (CL_req - CL0_est)/CL_alpha_rad;
        alpha_cmd = min(max(alpha_raw, alpha_min), alpha_max);
        al_sat_err = abs(alpha_raw - alpha_cmd);

        a_vert_req = dot(a_cmd_ecef_guess,u_up);
        a_h_req_vec = a_cmd_ecef_guess - a_vert_req*u_up;
        a_lat_req = dot(a_h_req_vec,right_h);

        phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
        phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

        theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
        theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);
        theta_cmd = theta_cmd(1);

        if R_to_T > 900e3
            M_cmd = M_cmd_far;
        elseif R_to_T > 250e3
            M_cmd = M_cmd_mid;
        else
            M_cmd = M_cmd_near;
        end
        V_cmd = max(M_cmd*a_snd, V_floor);

        [CL, CD, CT_now] = aero_coeffs(Ma, alpha_cmd, dT);
        L = CL*qbar*S_ref;
        D = CD*qbar*S_ref;
        CT_ff = D / max(qbar*S_ref,1);

        v_err = V_cmd - V;
        int_v_tmp = int_v + v_err*dt_use;
        int_v_tmp = min(max(int_v_tmp,-int_v_lim),int_v_lim);
        CT_cmd = CT_ff + Kpv*v_err + Kiv*int_v_tmp;
        CT_cmd = min(max(CT_cmd,CT_min),CT_max);

        dT_grid = linspace(dT_min,dT_max,31);
        CT_grid = zeros(size(dT_grid));
        for ii=1:numel(dT_grid)
            [~,~,CT_grid(ii)] = aero_coeffs(Ma, alpha_cmd, dT_grid(ii));
        end
        [~,ix] = min(abs(CT_grid - CT_cmd));
        dT_tmp = dT_grid(ix);

        [~,~,CT_tmp] = aero_coeffs(Ma, alpha_cmd, dT_tmp);
        T_eng = CT_tmp * qbar * S_ref;

        F_drag_e   = -D * (v/max(norm(v),1e-6));
        F_lift_e   =  L * u_up;
        if norm(v_h_vec) > 1e-6
            fwd = v_h_vec / norm(v_h_vec);
        else
            fwd = v / max(norm(v),1e-6);
        end
        F_thrust_e = T_eng * fwd;
        F_ecef = F_drag_e + F_lift_e + F_thrust_e;

        [lat_g, lon_g, h_g] = ecef2lla_cgcs2000(r');
        g_ecef = gravity_cgcs2000_ecef(lat_g, lon_g, h_g);

        a_ecef = g_ecef + F_ecef/m ...
             - 2*cross(omega_ie, v) ...
             - cross(omega_ie, cross(omega_ie, r));

        v_pred = v + a_ecef*dt_use;
        r_pred = r + v_pred*dt_use;

        [~, ~, h_pred] = ecef2lla_cgcs2000(r_pred');
        V_pred = max(norm(v_pred),1e-3);

        u_up_pred = r_pred / max(norm(r_pred),1e-6);
        v_up_pred = dot(v_pred, u_up_pred);
        v_h_pred = v_pred - v_up_pred*u_up_pred;
        Vh_pred = max(norm(v_h_pred),1e-6);
        gamma_pred = atan2(v_up_pred, Vh_pred);
        if ~isscalar(gamma_pred), gamma_pred = gamma_pred(1); end

        R_pred = norm(rT - r_pred);

        J = 1.0e-6*(R_pred^2) + ...
            8.0e-4*((h_pred - h_cmd_cruise)^2) + ...
            1.0e-4*((V_pred - V_cmd)^2) + ...
            1.0e-3*(a_lat_cmd^2) + ...
            5.0e-4*(gamma_pred^2) + ...
            2.0e-3*(al_sat_err^2);

        if J < bestJ
            bestJ = J;
            best_a_lat = a_lat_cmd;
            best_alpha  = alpha_cmd;
            best_phi    = phi_cmd;
            best_theta  = theta_cmd;
            best_CL     = CL;
            best_CD     = CD;
            best_CT     = CT_tmp;
            best_CTff   = CT_ff;
            best_CLreq  = CL_req;
            best_dT     = dT_tmp;
        end
    end

    %% ========== 使用最优横向控制 ==========
    a_lat_cmd = best_a_lat;
    alpha_cmd = best_alpha;
    phi_cmd   = best_phi;
    theta_cmd = best_theta;
    theta_cmd = theta_cmd(1);

    % 最终纵向CL需求再算一次，用于实际执行
    CL_req = best_CLreq;
    CL_req = min(max(CL_req, CL_floor), CL_ceil);

    % alpha饱和卸载
    alpha_sat = min(max(alpha_cmd, alpha_min), alpha_max);
    if abs(alpha_cmd - alpha_sat) > 1e-9
        over = abs(alpha_cmd - alpha_sat)/max(abs(alpha_max-alpha_min),1e-6);
        unload = min(1.0, alpha_unload_gain * (0.2 + 3.0*over));
        int_h = int_h * (1.0 - 0.5*unload);
    end
    alpha_cmd = alpha_sat;

    phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);
    theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);
    theta_cmd = theta_cmd(1);

    psi_cmd = atan2(dot(u_vh,eE), dot(u_vh,eN));
    if ~isscalar(psi_cmd), psi_cmd = psi_cmd(1); end

    phi = phi_cmd; theta = theta_cmd; psi = psi_cmd;

    %% ========== 速度控制 ==========
    if R_to_T > 900e3
        M_cmd = M_cmd_far;
    elseif R_to_T > 250e3
        M_cmd = M_cmd_mid;
    else
        M_cmd = M_cmd_near;
    end
    V_cmd = max(M_cmd*a_snd, V_floor);

    [CL, CD, CT_now] = aero_coeffs(Ma, alpha_cmd, dT);
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

    %% ========== 力与积分 ==========
    if norm(v_h_vec) > 1e-6
        fwd = v_h_vec / norm(v_h_vec);
    else
        fwd = v / max(norm(v),1e-6);
    end

    F_drag_e   = -D * (v/max(norm(v),1e-6));
    F_lift_e   =  L * u_up;
    F_thrust_e =  T_eng * fwd;
    F_ecef = F_drag_e + F_lift_e + F_thrust_e;

    [lat_g, lon_g, h_g] = ecef2lla_cgcs2000(r');
    g_ecef = gravity_cgcs2000_ecef(lat_g, lon_g, h_g);

    a_ecef = g_ecef + F_ecef/m ...
         - 2*cross(omega_ie, v) ...
         - cross(omega_ie, cross(omega_ie, r));

    v = v + a_ecef*dt_use;
    r = r + v*dt_use;
    t_now = t_now + dt_use;

    %% ========== 记录 ==========
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
    aLat_hist(k) = a_lat_cmd;
    CLreq_hist(k) = CL_req;
    alpha_hist(k) = alpha_cmd;
    phi_hist(k)   = phi_cmd;
    theta_hist(k) = theta_cmd;
    T_hist(k)     = T_eng;
    D_hist(k)     = D;
    L_hist(k)     = L;

    if any(~isfinite([r;v;Ma;qbar;CL;CD;CT;dT]))
        stop_reason = "NaN/Inf";
        fprintf('NaN/Inf at t=%.2f s\n', t_now);
        break;
    end

    if mod(round(t_now,1),10.0)==0
        fprintf("t=%.0f  R=%.0fkm  h=%.0f  V=%.0f  a_lat=%.1f  CLreq=%.2f  alpha=%.2fdeg\n", ...
            t_now, R_to_T/1000, h, V, a_lat_cmd, CL_req, alpha_cmd*180/pi);
    end

    k = k + 1;
end

%% ========== 截断有效数据 ==========
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
aLatk = aLat_hist(valid);
CLreqk = CLreq_hist(valid);
alphak = alpha_hist(valid);
phik   = phi_hist(valid);
thetak = theta_hist(valid);
Tk     = T_hist(valid);
Dforce = D_hist(valid);
Lforce = L_hist(valid);

fprintf('Simulation stop reason: %s, t=%.2f s, min range=%.1f m\n', stop_reason, tt(end), min(Dk));

%% ========== 绘图 ==========
figure;
plot(tt, Dk/1000, 'm', 'LineWidth',1.7); grid on;
xlabel('Time (s)'); ylabel('Distance to Target (km)');
title('Distance to Target');

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
plot(tt, Hk/1000, 'LineWidth',1.4); grid on; yline(h_cmd_cruise/1000,'r--');
xlabel('Time (s)'); ylabel('Altitude (km)');
title('Altitude Hold');

figure;
plot(tt, CLreqk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('CL_{req}');
title('Required Lift Coefficient');

figure;
plot(tt, alphak*180/pi, 'LineWidth',1.3); grid on;
yline(alpha_min*180/pi,'r--'); yline(alpha_max*180/pi,'r--');
xlabel('Time (s)'); ylabel('\alpha (deg)');
title('Angle of Attack (limited)');

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
plot(tt, CTk, 'k', 'LineWidth',1.2); grid on;
xlabel('Time (s)'); ylabel('C_T (actual)');
title('CT Actual');

figure;
plot(tt, CLk, 'b', 'LineWidth',1.2); hold on; grid on;
plot(tt, CDk, 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Coefficient');
title('C_L / C_D');
legend('C_L','C_D','Location','best');

figure;
plot(tt, aLatk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('a_{lat,cmd} (m/s^2)');
title('Lateral Acceleration Command');

figure;
plot(tt, Tk, 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('Thrust (N)');
title('Engine Thrust (T = CT*q*S)');

figure;
plot(tt, (Tk - Dforce), 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('(T - D) (N)');
title('Excess Thrust');

figure;
plot3(Rk(:,1)/1e3, Rk(:,2)/1e3, Rk(:,3)/1e3, 'b','LineWidth',1.3); hold on; grid on; axis equal;
[xe,ye,ze] = sphere(60);
surf(Re*xe/1e3, Re*ye/1e3, Re*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none');
plot3(rT(1)/1e3, rT(2)/1e3, rT(3)/1e3, 'ro','MarkerFaceColor','r');
xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
title('Trajectory in ECEF');
legend('Vehicle','Earth','Target','Location','best');

%% ========== 发射系 ENU ==========
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
title('Trajectory in Launch Frame (ENU)');

%% ========== 轨迹经纬度图 ==========
Ntraj = size(Rk,1);
lat_traj = nan(Ntraj,1);
lon_traj = nan(Ntraj,1);
h_traj   = nan(Ntraj,1);

for ii = 1:Ntraj
    [lat_traj(ii), lon_traj(ii), h_traj(ii)] = ecef2lla_cgcs2000(Rk(ii,:));
end

figure;
plot(lon_traj, lat_traj, 'b', 'LineWidth', 1.5); grid on; hold on;
plot(lon0, lat0, 'go', 'MarkerFaceColor','g');
plot(lonT, latT, 'ro', 'MarkerFaceColor','r');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
title('Trajectory Ground Track (Lat/Lon)');
legend('Trajectory','Launch','Target','Location','best');

figure;
plot3(lon_traj, lat_traj, h_traj/1000, 'b', 'LineWidth', 1.5); grid on; hold on;
plot3(lon0, lat0, h0/1000, 'go', 'MarkerFaceColor','g');
plot3(lonT, latT, hT/1000, 'ro', 'MarkerFaceColor','r');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); zlabel('Altitude (km)');
title('Trajectory in LLA (Lon-Lat-Alt)');
legend('Trajectory','Launch','Target','Location','best');
view(3);

figure;
subplot(2,1,1);
plot(lon_traj, h_traj/1000, 'LineWidth', 1.3); grid on;
xlabel('Longitude (deg)'); ylabel('Altitude (km)'); title('Altitude vs Longitude');

subplot(2,1,2);
plot(lat_traj, h_traj/1000, 'LineWidth', 1.3); grid on;
xlabel('Latitude (deg)'); ylabel('Altitude (km)'); title('Altitude vs Latitude');

%% ====== local function ======
function CL = getCL(Ma, alpha, dT)
    [CL,~,~] = aero_coeffs(Ma, alpha, dT);
end