clear;
clc;

%% ===================== 时间 =====================
dt = 0.1;
T  = 2000;
t  = 0:dt:T;
N  = length(t);

%% ===================== 地球常量 =====================
mu  = 3.986004418e14;
Re  = 6378137.0;

% Debug物理开关（C）
use_simple_gravity = true;
use_rotation_terms = false;

%% ===================== 飞行器参数 =====================
m     = 671.33;
S_ref = 0.2986;

%% ===================== 目标 =====================
h_target = 30e3;
lat_T = 13.35;
lon_T = 144.55;
rT = lla2ecef_cgcs2000(lat_T, lon_T, h_target); rT = rT(:);
vT = [0;0;0];

%% ===================== 初始条件 =====================
h_init = 30e3;
lat0 = 19.2;
lon0 = 110.5;
r = lla2ecef_cgcs2000(lat0, lon0, h_init); r = r(:);

Ma_init = 6.5;
[a0,~] = atmos_simple(h_init);
V0 = Ma_init * a0;

u_up0 = r / norm(r);
v_dir = rT - (rT' * u_up0) * u_up0;
v_dir = v_dir / max(norm(v_dir),1e-6);
v = V0 * v_dir;

% 姿态初始化（用于alpha/beta估算）
eU = [cosd(lat0)*cosd(lon0); cosd(lat0)*sind(lon0); sind(lat0)];
eN = [-sind(lat0)*cosd(lon0); -sind(lat0)*sind(lon0); cosd(lat0)];
x_h = (rT-r) - ((rT-r)'*(r/norm(r))) * (r/norm(r));
if norm(x_h) < 1e-6, xb = eN; else, xb = x_h / norm(x_h); end
yb = cross(eU, xb); yb = yb / max(norm(yb),1e-6);
zb = cross(xb, yb); zb = zb / max(norm(zb),1e-6);
C_b2e_0 = [xb, yb, zb];

theta = asin(-C_b2e_0(3,1));
psi   = atan2(C_b2e_0(2,1), C_b2e_0(1,1));
phi   = atan2(C_b2e_0(3,2), C_b2e_0(3,3));

% 舵面/油门
de1 = 0; de2 = 0; dr = 0; dT = 0.3;

%% ===================== 控制参数 =====================
Npn = 2.8;
a_cmd_max = 80;

k_psi = 1.8;
psi_rate_max = 8*pi/180;

phi_lim   = 60*pi/180;
theta_lim = 30*pi/180;
alpha_min = -2*pi/180;
alpha_max = 15*pi/180;

Kph = 0.0010;
Kdh = 0.06;

% 速度指令
V_floor  = 700;
tau_vcmd = 15;
V_cmd = V0;

% -------- 推力优先：CT前馈 + PI速度环 --------
Kp_v = 0.0018;          % PI比例
Ki_v = 0.00025;         % PI积分
int_v_err = 0;          % 积分状态
int_v_lim = 8000;       % 防积分饱和

dT_min = 0.18;          % 提高最小油门
dT_max = 1.00;

%% ===================== 记录 =====================
Rhist = nan(N,3);
Vhist = nan(N,3);
Hhist = nan(N,1);
Dist_hist = nan(N,1);
Vc_hist = nan(N,1);
Vmag_hist = nan(N,1);
Vcmd_hist = nan(N,1);

dT_hist = nan(N,1);
CTff_hist = nan(N,1);
CTcmd_hist = nan(N,1);

CL_hist = nan(N,1);
CL_req_hist = nan(N,1);
CD_hist = nan(N,1);

alpha_hist = nan(N,1);
alpha_cmd_hist = nan(N,1);

a_up_aero_hist  = nan(N,1);
a_up_total_hist = nan(N,1);

stop_k = N;
stop_reason = "completed";

%% ===================== 主循环 =====================
for k = 1:N
    % 地理量
    [lat, lon, h] = ecef2lla_cgcs2000(r');
    Hhist(k) = h;

    if h <= 25e3
        stop_k = k; stop_reason = "terrain impact";
        fprintf('Vehicle impacted terrain at t=%.1f s\n', t(k));
        break;
    end

    % 大气
    v_air = v;
    Vair = max(norm(v_air),1e-3);
    [a, rho] = atmos_simple(max(h,0));
    Ma = max(min(Vair/max(a,1e-3), 7.0), 0.0);
    qbar = 0.5*rho*Vair^2;

    % 姿态矩阵（仅用于alpha/beta估算）
    cph=cos(phi); sph=sin(phi); cth=cos(theta); sth=sin(theta); cps=cos(psi); sps=sin(psi);
    Cned_b = [cth*cps, sph*sth*cps-cph*sps, cph*sth*cps+sph*sps;
              cth*sps, sph*sth*sps+cph*cps, cph*sth*sps-sph*cps;
              -sth,    sph*cth,             cph*cth];

    slat = sin(lat); clat = cos(lat);
    slon = sin(lon); clon = cos(lon);
    eN_e = [-slat*clon; -slat*slon; clat];
    eE_e = [-slon;      clon;       0];
    eD_e = [-clat*clon; -clat*slon; -slat];
    Ce_ned = [eN_e, eE_e, eD_e];
    Ce_b = Ce_ned * Cned_b;
    Cb_e = Ce_b';

    % alpha,beta估算
    v_air_b = Cb_e * v_air;
    u = v_air_b(1); vv = v_air_b(2); w = v_air_b(3);
    alpha = atan2(w,u);
    beta  = asin(max(min(vv/Vair,1),-1));
    alpha = max(min(alpha, 8*pi/180), -8*pi/180);
    beta  = max(min(beta,  5*pi/180), -5*pi/180);

    % 重力（C）
    if use_simple_gravity
        g_ecef = -mu/norm(r)^3 * r;
    end

    %% ---------- 导引 ----------
    r_rel = rT - r;
    R_dist = max(norm(r_rel),1);
    u_los = r_rel / R_dist;
    v_rel = vT - v;
    Vc = -dot(r_rel, v_rel) / R_dist;
    omega_los = cross(r_rel, v_rel) / (R_dist^2 + 1e-6);

    if R_dist < 100e3, Npn_use = 1.5;
    elseif R_dist < 500e3, Npn_use = 2.0;
    else, Npn_use = Npn;
    end

    a_pn_ecef = Npn_use * max(Vc,0) * cross(u_los, omega_los);
    if norm(a_pn_ecef) > a_cmd_max
        a_pn_ecef = a_pn_ecef / norm(a_pn_ecef) * a_cmd_max;
    end

    if use_rotation_terms
        omega_ie = [0;0;7.2921150e-5];
        a_req_ecef = a_pn_ecef - g_ecef + 2*cross(omega_ie,v) + cross(omega_ie,cross(omega_ie,r));
    else
        a_req_ecef = a_pn_ecef - g_ecef;
    end

    % 局部基底
    u_up = r / norm(r);
    eE_now = [-sin(lon); cos(lon); 0];
    eN_now = [-sin(lat)*cos(lon); -sin(lat)*sin(lon); cos(lat)];

    % 水平速度方向
    v_h = v_air - dot(v_air,u_up)*u_up;
    if norm(v_h) > 1e-6
        v_h_hat = v_h / norm(v_h);
        chi_v = atan2(dot(v_h_hat,eE_now), dot(v_h_hat,eN_now));
    else
        v_h_hat = eN_now;
        chi_v = psi;
    end

    % LOS航向
    los_h = u_los - dot(u_los,u_up)*u_up;
    if norm(los_h) > 1e-8
        los_h_hat = los_h / norm(los_h);
        chi_los = atan2(dot(los_h_hat,eE_now), dot(los_h_hat,eN_now));
    else
        chi_los = chi_v;
    end

    % 航向控制
    e_psi = wrapToPi(chi_los - psi);
    psi_rate_cmd = max(min(k_psi*e_psi, psi_rate_max), -psi_rate_max);
    psi_cmd = wrapToPi(psi + psi_rate_cmd*dt);

    % 垂直/水平分配
    a_vert_geo = dot(a_req_ecef, u_up);
    a_h_geo = a_req_ecef - a_vert_geo*u_up;
    a_h = norm(a_h_geo);

    dh = h_init - h;
    vh = dot(v_air, u_up);
    a_hold = Kph*dh - Kdh*vh;
    a_vert_cmd = a_vert_geo + a_hold;

    CL_max = 1.2;
    a_L_cap = (CL_max*qbar*S_ref)/m;

    if R_dist < 200e3
        vert_share = 0.45; 
    else
        vert_share = 0.65;
    end

    a_vert_cmd = max(min(a_vert_cmd, vert_share*a_L_cap), -vert_share*a_L_cap);

    a_h_cap = sqrt(max(a_L_cap^2 - a_vert_cmd^2,0));
    if a_h > a_h_cap && a_h > 1e-6
        a_h_geo = a_h_geo * (a_h_cap/a_h);
        a_h = a_h_cap;
    end

    % 滚转
    right_dir = cross(u_up, v_h_hat); right_dir = right_dir / max(norm(right_dir),1e-6);
    a_lat = dot(a_h_geo, right_dir);
    phi_cmd = atan2(a_lat, max(abs(a_vert_cmd),1.0));
    phi_cmd = max(min(phi_cmd, phi_lim), -phi_lim);

    % alpha_cmd
    a_L_req = sqrt(a_vert_cmd^2 + a_h^2);
    CL_req = (m*a_L_req) / max(qbar*S_ref,1);

    de1_deg = de1*180/pi; de2_deg = de2*180/pi;
    CL0 = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.006515*(de1_deg + de2_deg);
    CL_alpha_rad = (0.07235 - 0.003368*Ma) * (180/pi);
    CL_alpha_rad = max(CL_alpha_rad, 0.1);

    alpha_cmd = (CL_req - CL0) / CL_alpha_rad;
    alpha_cmd = max(min(alpha_cmd, alpha_max), alpha_min);

    % 俯仰
    v_vert = dot(v_air, u_up);
    v_horz = norm(v_air - v_vert*u_up);
    gamma_now = atan2(v_vert, v_horz);
    theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
    theta_cmd = max(min(theta_cmd, theta_lim), -theta_lim);

    % 姿态更新
    phi = phi_cmd; theta = theta_cmd; psi = psi_cmd;

    %% ---------- 速度指令 ----------
    if R_dist > 900e3
        Vcmd_base = 6.0*a;
    elseif R_dist > 400e3
        Vcmd_base = 4.8*a;
    elseif R_dist > 150e3
        Vcmd_base = 3.8*a;
    else
        Vcmd_base = 3.2*a;
    end
    Vcmd_base = max(Vcmd_base, V_floor);
    V_cmd = V_cmd + (Vcmd_base - V_cmd)*dt/tau_vcmd;

    %% ---------- 气动 ----------
    [CL,CD,CY,~,~,~] = aero_coeffs(Ma,alpha_cmd,beta,de1,de2,dr,dT);
    CY = 0;   % B: 关侧力

    L = CL*qbar*S_ref;
    D = CD*qbar*S_ref;
    CT_now = thrust_coeffs(Ma,alpha_cmd,beta,dT);

    % 前向方向（ECEF）
    if norm(v_h) > 1e-6
        fwd = v_h / norm(v_h);
    else
        fwd = v_air / max(norm(v_air),1e-6);
    end
    right = cross(u_up,fwd); right = right/max(norm(right),1e-6);

    %% ---------- 推力优先控制（CT前馈 + PI） ----------
    % 1) 前馈：用阻力平衡估计所需CT
    %    T_req ~= D + m*a_x_des, 其中 a_x_des 用速度误差生成
    a_x_des = 0.02*(V_cmd - Vair);                        % 期望沿速度加速度
    T_req_ff = D + m*a_x_des;
    CT_ff = T_req_ff / max(qbar*S_ref, 1.0);

    % 2) PI纠偏（对CT增量）
    v_err = V_cmd - Vair;
    int_v_err = int_v_err + v_err*dt;
    int_v_err = max(min(int_v_err, int_v_lim), -int_v_lim);

    dCT_pi = Kp_v*v_err + Ki_v*int_v_err;
    CT_cmd = CT_ff + dCT_pi;

    % 3) 由CT_cmd反求dT（单调搜索，避免未知反函数）
    dT_grid = linspace(dT_min, dT_max, 41);
    CT_grid = zeros(size(dT_grid));
    for ii = 1:numel(dT_grid)
        CT_grid(ii) = thrust_coeffs(Ma,alpha_cmd,beta,dT_grid(ii));
    end
    [~,ii_best] = min(abs(CT_grid - CT_cmd));
    dT_unsat = dT_grid(ii_best);

    % 4) 油门一阶执行器
    tau_dT = 0.6;
    dT = dT + (dT_unsat - dT)*dt/tau_dT;
    dT = max(min(dT, dT_max), dT_min);

    % 5) anti-windup: 若饱和且误差继续推向饱和，回退积分
    at_upper = (dT >= dT_max-1e-6) && (v_err > 0);
    at_lower = (dT <= dT_min+1e-6) && (v_err < 0);
    if at_upper || at_lower
        int_v_err = int_v_err - v_err*dt*0.7;
    end

    % 实际推力
    CT = thrust_coeffs(Ma,alpha_cmd,beta,dT);
    T_eng = CT*qbar*S_ref;

    %% ---------- A: ECEF直接构力 ----------
    F_drag_e   = -D * (v_air / max(norm(v_air),1e-6));
    F_lift_e   =  L * u_up;
    F_side_e   =  0 * right; 
    F_thrust_e =  T_eng * fwd;

    F_ecef = F_drag_e + F_lift_e + F_thrust_e;

    a_aero_ecef = (F_drag_e + F_lift_e)/m;
    a_up_aero = dot(a_aero_ecef, u_up);

    a_total = F_ecef/m + g_ecef;
    a_up_total = dot(a_total, u_up);

    %% ---------- 平动 ----------
    if use_rotation_terms
        omega_ie = [0;0;7.2921150e-5];
        a_ecef = g_ecef + F_ecef/m - 2*cross(omega_ie,v) - cross(omega_ie,cross(omega_ie,r));
    else
        a_ecef = g_ecef + F_ecef/m;
    end

    v = v + a_ecef*dt;
    r = r + v*dt;

    %% ---------- 记录 ----------
    Rhist(k,:) = r.';
    Vhist(k,:) = v.';
    Dist_hist(k) = R_dist;
    Vc_hist(k) = Vc;
    Vmag_hist(k) = Vair;
    Vcmd_hist(k) = V_cmd;

    dT_hist(k) = dT;
    CTff_hist(k) = CT_ff;
    CTcmd_hist(k) = CT_cmd;

    CL_hist(k) = CL;
    CL_req_hist(k) = CL_req;
    CD_hist(k) = CD;

    alpha_hist(k) = alpha;
    alpha_cmd_hist(k) = alpha_cmd;

    a_up_aero_hist(k) = a_up_aero;
    a_up_total_hist(k) = a_up_total;

    if R_dist < 5e3
        stop_k = k; stop_reason = "target reached";
        fprintf('Target reached at t=%.1f s, miss distance=%.1f m\n', t(k), R_dist);
        break;
    end

    if any(~isfinite([r;v;phi;theta;psi;Ma;qbar;alpha;beta;CL;CD;CT]))
        stop_k = k; stop_reason = "NaN/Inf";
        fprintf('NaN/Inf at step %d, t=%.1f s\n', k, t(k));
        break;
    end
end

%% ===================== 有效数据 =====================
valid = isfinite(Dist_hist) & (Dist_hist > 0);
tt = t(valid);

fprintf('Simulation stop reason: %s, t=%.1f s\n', stop_reason, t(stop_k));

%% ===================== 绘图 =====================
figure;
plot(tt, Dist_hist(valid)/1000, 'm-', 'LineWidth',1.6); grid on;
xlabel('Time (s)'); ylabel('Distance to Target (km)');
title('Distance to Target (CT-FF + PI)');

figure;
plot(tt, Vc_hist(valid), 'b-', 'LineWidth',1.4); grid on; yline(0,'r--');
xlabel('Time (s)'); ylabel('V_c (m/s)');
title('Closing Speed (CT-FF + PI)');

figure;
plot(tt, Vmag_hist(valid), 'k-', 'LineWidth',1.4); hold on; grid on;
plot(tt, Vcmd_hist(valid), 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (m/s)');
title('Speed Tracking (CT-FF + PI)');
legend('|V|','V_{cmd}','Location','best');

figure;
plot(tt, dT_hist(valid), 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('dT');
title('Throttle Command');

figure;
plot(tt, CTff_hist(valid), 'b--', 'LineWidth',1.2); hold on; grid on;
plot(tt, CTcmd_hist(valid), 'r-', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('C_T');
title('CT Feedforward / CT Command');
legend('CT_{ff}','CT_{cmd}','Location','best');

figure;
plot(tt, Hhist(valid)/1000, 'LineWidth',1.4); grid on;
xlabel('Time (s)'); ylabel('Altitude (km)');
title('Altitude');

figure;
plot(tt, a_up_aero_hist(valid), 'b-', 'LineWidth',1.3); hold on; grid on;
plot(tt, a_up_total_hist(valid), 'm-', 'LineWidth',1.3);
yline(0,'k--');
xlabel('Time (s)');
ylabel('Acceleration along +Up (m/s^2)');
title('Up-axis Acceleration Diagnostics');
legend('a_{up,aero}','a_{up,total}','0','Location','best');

figure;
plot3(Rhist(valid,1)/1e3, Rhist(valid,2)/1e3, Rhist(valid,3)/1e3, 'b', 'LineWidth',1.2); hold on; grid on; axis equal;
[xe,ye,ze] = sphere(60);
surf(Re*xe/1e3, Re*ye/1e3, Re*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none','FaceColor',[0.2 0.6 1.0]);
plot3(rT(1)/1e3, rT(2)/1e3, rT(3)/1e3, 'ro', 'MarkerFaceColor','r');
xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
title('Trajectory in ECEF (CT-FF + PI)');
legend('Vehicle','Earth','Target','Location','best');