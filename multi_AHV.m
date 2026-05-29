% 4机编队 3DOF（ECEF）仿真：基于单机 multi_AHV 扩展
% - 领机(1)追真实目标 rT
% - 僚机(2/3/4)追"虚拟结构(VS)"队形目标点 rT_i = r1 + ex*x_i + ey*y_i
% - 队形为V型/菱形拓扑：1-2,1-3,2-4,3-4 为主要邻接
% - 收窄策略：以领机到目标距离 R 为自变量，从远端 d_adj=3000m 逐步收窄到近端 d_adj=1500m
% - 约束：横向展开 <= 2km (|y_wing|<=1000m), 最近邻安全距离 >= 1km（用排斥加速度屏障）
%
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
nVeh = 4;
m = 671.33 * ones(nVeh,1);
S_ref = 0.2986 * ones(nVeh,1);

%% ================= 起点/终点设置（沿用单机） =================
h0 = 30e3; lat0 = 19.2; lon0 = 110.5;   % lat/lon 单位：度
r0 = lla2ecef_cgcs2000(lat0, lon0, h0); r0 = r0(:);

latT = 13.35; lonT = 144.55; hT = 30e3; % 度
rT = lla2ecef_cgcs2000(latT, lonT, hT); rT = rT(:);

[a0,~] = atmos_simple(h0);
V0 = 6.5 * a0;

%% ================= 航点（沿用单机：两点航段） =================
wps_lla = [lat0, lon0, h0;
    latT, lonT, hT];
nWP = size(wps_lla,1);
wps_ecef = zeros(nWP,3);
for i=1:nWP
    wps_ecef(i,:) = lla2ecef_cgcs2000(wps_lla(i,1), wps_lla(i,2), wps_lla(i,3));
end
wp_switch_dist = 100e3;

%% ================= 编队参数：收窄策略（用户给定） =================
% 收窄起始距离（领机到目标）
R_narrow_start = 1000e3;   % 1000 km

% 近端邻接目标距离（用户给定）
d_adj_near = 1500;         % m

% 远端邻接目标距离（初始队形）
d_adj_far  = 3000;         % m

% 横向展开限制：左右翼机 |y|<=1000 m  -> 展开 <= 2 km
wing_y_max = 1000;         % m

% 最近邻安全距离（硬安全）
d_safe_min = 1000;         % m

% 近端收窄结束距离（为了平滑，这里取 120 km；你也可按需求调整）
R_narrow_end = 120e3;      % m

% 收窄平滑（0~1）函数：R从 start->end 变为 s=0->1
narrow_smooth = @(R) min(max((R_narrow_start - R) / max(R_narrow_start - R_narrow_end, 1), 0), 1);

%% ================= 编队控制（严格队形） =================
% 说明：
% - 领机：沿原始航点航段 L1 追 rT
% - 僚机：不再沿航段 L1，而是对队形目标点 rT_i 做 L1/LOS 横向制导（强闭环），以严格保持菱形/V
% - 为避免和安全屏障冲突，僚机横向指令限幅更严格一些

K_form_xt = 1.0;           % 僚机横向误差增益（可调 0.5~2.0）
a_lat_max_wing = 90;       % 僚机横向加速度限幅（m/s^2）

%% ================= 初始编队布阵（V型/菱形拓扑，邻接约3km） =================
% 用领机局部水平基：ex 指向目标水平前向，ey 水平侧向（左），ez 天顶
u_up0 = r0/norm(r0);
ex0 = rT - r0;
ex0 = ex0 - dot(ex0,u_up0)*u_up0;
ex0 = ex0 / max(norm(ex0),1e-6);
ey0 = cross(u_up0, ex0);  ey0 = ey0 / max(norm(ey0),1e-6);
ez0 = u_up0;

% 初始邻接距离采用远端 d_adj_far：等边三角高度 h = sqrt(3)/2*d
h_far = sqrt(3)/2 * d_adj_far;

% 1: 头机；2/3: 左/右翼（距头机d）；4: 后机（距左右翼d）
offV0 = [ ...
    0,        0,      0;          % 1
   -h_far,   -d_adj_far/2, 0;     % 2
   -h_far,    d_adj_far/2, 0;     % 3
   -2*h_far,  0,      0];         % 4

r = zeros(3,nVeh);
for i=1:nVeh
    r(:,i) = r0 + offV0(i,1)*ex0 + offV0(i,2)*ey0 + offV0(i,3)*ez0;
end

%% ================= 初始速度方向设置（沿用单机：朝目标水平分量） =================
seg0 = wps_ecef(2,:).' - wps_ecef(1,:).';
seg0_h = seg0 - dot(seg0,u_up0)*u_up0;
if norm(seg0_h) < 1e-6
    seg0_h = seg0;
end
u_vh0 = seg0_h / norm(seg0_h);

v = zeros(3,nVeh);
for i=1:nVeh
    v(:,i) = V0 * u_vh0;
end

%% ================= 控制参数设置（沿用单机） =================
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

int_h = zeros(nVeh,1);

phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2;
alpha_max =  10;

alpha_unload_gain = 0;   % 如需饱和卸载，可改为 0.5~2.0

% 速度控制
Kpv = 0.0018;
Kiv = 0.00025;
int_v = zeros(nVeh,1);
int_v_lim = 8e3;

CT_min = 0.025;
CT_max = 0.3;
dT_min = 0.01;
dT_max = 1.00;
tau_dT = 0.3;
dT = 0.5 * ones(nVeh,1);

%% ================= 固定网格（加速：循环内复用，不重复分配） =================
dT_grid = linspace(dT_min, dT_max, 41);

%% ================= 初始配平（按单机方法，统一给所有机） =================
phi = zeros(nVeh,1); theta = zeros(nVeh,1); psi = zeros(nVeh,1);

use_trim_init = true;
if use_trim_init
    [a0t, rho0] = atmos_simple(h0);
    Ma0   = V0 / max(a0t,1e-3);
    qbar0 = 0.5 * rho0 * V0^2;

    R0 = 6371000;
    g = g0 * (R0/(R0+h0))^2;
    CL_req0 = m(1)*g / (qbar0*S_ref(1));

    alpha0 = fminbnd(@(a) (aero_coeffs(Ma0, a, dT(1)) - CL_req0)^2, alpha_min, alpha_max);

    [~, CD0, ~] = aero_coeffs(Ma0, alpha0, dT(1));
    D0 = CD0 * qbar0 * S_ref(1);
    CT_req0 = D0 / max(qbar0*S_ref(1),1);

    CT_grid0 = zeros(size(dT_grid));
    for ii=1:numel(dT_grid)
        [~,~,CT_grid0(ii)] = aero_coeffs(Ma0, alpha0, dT_grid(ii));
    end
    [~,ix] = min(abs(CT_grid0 - CT_req0));
    dT(:) = dT_grid(ix);

    eE0 = [-sind(lon0); cosd(lon0); 0];
    eN0 = [-sind(lat0)*cosd(lon0); -sind(lat0)*sind(lon0); cosd(lat0)];
    chi0 = atan2(dot(u_vh0,eE0), dot(u_vh0,eN0));

    gamma0 = 0;
    phi(:)   = 0;
    theta(:) = gamma0 + deg2rad(alpha0);
    psi(:)   = chi0;

    fprintf('Trim init (formation): alpha0=%.3f deg, dT0=%.3f, theta0=%.3f deg\n', ...
        alpha0, dT(1), theta(1)*180/pi);
end

%% ================= 马赫指令（沿用单机） =================
M_cmd_far  = 6.5;
M_cmd_mid  = 5.5;
M_cmd_near = 5.0;
V_floor    = 1100;

%% ================= 终端捕获参数（以领机为准） =================
R_hit            = 1000;
R_term1          = 120e3;
R_term2          = 40e3;
R_pass_check     = 20e3;
R_abort_diverge  = 50e3;

%% ================= 僚机航点索引（保持与单机一致的L1框架；领机使用） =================
wp_idx = 2 * ones(nVeh,1);
r_wp_prev = repmat(wps_ecef(1,:).', 1, nVeh);
r_wp_next = repmat(wps_ecef(2,:).', 1, nVeh);

%% ================= 队形安全屏障参数 =================
% 屏障只在 d<d_safe_min 时介入；a_rep_max 限幅，避免强烈震荡
barrier_on = true;
a_rep_max = 12;        % m/s^2
k_rep = 1.5;           % 排斥强度（配合a_rep_max使用）

%% ================= 记录数组（4机维度） =================
t_log = nan(N_max,1);
Rhist = nan(N_max,3,nVeh);
Vhist = nan(N_max,3,nVeh);
Hhist = nan(N_max,nVeh);
Dist_hist = nan(N_max,nVeh);

Vmag_hist = nan(N_max,nVeh);
Vcmd_hist = nan(N_max,nVeh);
dT_hist   = nan(N_max,nVeh);
alpha_hist = nan(N_max,nVeh);

% 编队关键距离日志（用于验证约束）
d12_hist = nan(N_max,1);
d13_hist = nan(N_max,1);
d24_hist = nan(N_max,1);
d34_hist = nan(N_max,1);
width_hist = nan(N_max,1);

stop_reason = "completed";
k = 1;
t_now = 0;

prev_R = inf;
min_R  = inf;

%% ================= 可选：进度输出 =================
show_progress = true;
progress_every_k = 2000;

%% ================= 主循环 =================
while k <= N_max && t_now <= T_end

    % -------- 领机到目标距离（用于步长与收窄调度） --------
    r1 = r(:,1);
    v1 = v(:,1);
    r_rel_T1 = rT - r1;
    R_to_T1  = norm(r_rel_T1);
    min_R = min(min_R, R_to_T1);

    if show_progress && mod(k, progress_every_k)==0
        [~,~,h1_dbg] = ecef2lla_cgcs2000(r1');
        fprintf('[progress] k=%d, t=%.2f s, R_leader=%.0f km, h_leader=%.0f m\n', k, t_now, R_to_T1/1e3, h1_dbg);
        drawnow limitrate;
    end

    % 终止条件以领机为准
    [~,~,h1] = ecef2lla_cgcs2000(r1');
    if h1 <= 20e3
        stop_reason = "terrain impact (leader)";
        fprintf('Leader impacted terrain at t=%.1f s\n', t_now);
        break;
    end

    u_los_T1 = r_rel_T1 / max(R_to_T1,1);
    Vc_toward1 = dot(v1, u_los_T1);

    if R_to_T1 <= R_hit
        stop_reason = "target reached (leader)";
        fprintf('Leader reached target at t=%.1f s, miss=%.1f m\n', t_now, R_to_T1);
        break;
    end

    if (R_to_T1 < R_pass_check) && (Vc_toward1 <= 0)
        stop_reason = "passed closest approach (leader)";
        fprintf('Leader passed closest approach at t=%.1f s, min range=%.1f m\n', t_now, min_R);
        break;
    end

    if (prev_R < R_abort_diverge) && (R_to_T1 > prev_R + 30)
        stop_reason = "diverging after near pass (leader)";
        fprintf('Leader diverging after near pass at t=%.1f s, range=%.1f m\n', t_now, R_to_T1);
        break;
    end
    prev_R = R_to_T1;

    % 步长（按领机距离）
    if R_to_T1 < R_term2
        dt_use = 0.001;
    elseif R_to_T1 < R_term1
        dt_use = 0.01;
    else
        dt_use = dt_base;
    end

    % -------- 形成"队形参考坐标系"：ex 指向目标水平前向，ey 水平侧向 --------
    u_up = r1 / norm(r1);
    ex = rT - r1;
    ex = ex - dot(ex,u_up)*u_up;
    ex = ex / max(norm(ex),1e-6);
    ey = cross(u_up, ex);
    ey = ey / max(norm(ey),1e-6);

    % -------- 收窄调度：邻接距离从 3000 -> 1500，且横向 |y|<=1000 --------
    s = narrow_smooth(R_to_T1);    % 0(远)->1(近)
    d_adj = (1-s)*d_adj_far + s*d_adj_near;
    d_adj = max(d_adj, d_safe_min);

    y_w = min(0.5*d_adj, wing_y_max);
    x_w = sqrt(max(d_adj^2 - y_w^2, 0));
    x4 = -x_w - sqrt(max(d_adj^2 - y_w^2, 0));

    p_des = zeros(3,nVeh);
    p_des(:,1) = [0;0;0];
    p_des(:,2) = (-x_w)*ex + (-y_w)*ey;
    p_des(:,3) = (-x_w)*ex + ( y_w)*ey;
    p_des(:,4) = ( x4 )*ex;

    rT_veh = zeros(3,nVeh);
    rT_veh(:,1) = rT;
    for iv=2:nVeh
        rT_veh(:,iv) = r1 + p_des(:,iv);
    end

    % -------- 每机控制与积分 --------
    for iv = 1:nVeh
        ri = r(:,iv);
        vi = v(:,iv);
        rTi = rT_veh(:,iv);

        [lat_deg, lon_deg, h] = ecef2lla_cgcs2000(ri');

        if h <= 20e3
            stop_reason = sprintf('terrain impact (veh %d)', iv);
            fprintf('Vehicle %d impacted terrain at t=%.1f s\n', iv, t_now);
            k = N_max + 1;
            break;
        end

        % 大气
        V = max(norm(vi),1e-3);
        [a_snd, rho] = atmos_simple(max(h,0));
        Ma = max(min(V/max(a_snd,1e-3),8.0),0.0);
        qbar = 0.5*rho*V^2;

        % 局部基
        u_up_i = ri / norm(ri);
        eE = [-sind(lon_deg); cosd(lon_deg); 0];
        eN = [-sind(lat_deg)*cosd(lon_deg); -sind(lat_deg)*sind(lon_deg); cosd(lat_deg)];

        v_h = vi - dot(vi,u_up_i)*u_up_i;
        Vh = max(norm(v_h),1e-6);
        u_vh = v_h / Vh;

        % 航点切换（只对领机/保持兼容；僚机横向不使用航段）
        if norm(r_wp_next(:,iv) - ri) < wp_switch_dist && wp_idx(iv) < nWP
            wp_idx(iv) = wp_idx(iv) + 1;
            r_wp_prev(:,iv) = r_wp_next(:,iv);
            r_wp_next(:,iv) = wps_ecef(wp_idx(iv),:).';
        end

        %% ============= 横向制导 =============
        % 领机：沿航段 L1
        % 僚机：对队形目标点 rTi 做 L1/LOS（严格队形）

        if iv == 1
            seg = r_wp_next(:,iv) - r_wp_prev(:,iv);
            t_hat = seg / max(norm(seg),1e-6);
            t_h = t_hat - dot(t_hat,u_up_i)*u_up_i;
            t_h = t_h / max(norm(t_h),1e-6);

            r_from_start = ri - r_wp_prev(:,iv);
            r_xt = r_from_start - dot(r_from_start,t_h)*t_h;
            xt = norm(r_xt);

            if xt > 1e-6
                n_xt = r_xt/xt;
            else
                n_xt = cross(u_up_i,t_h);
                n_xt = n_xt/max(norm(n_xt),1e-6);
            end

            r_rel_T = rTi - ri;
            R_to_T = norm(r_rel_T);

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

            r_L1 = ri + L1_dist*t_h - min(xt,0.5*L1_dist)*n_xt;
            u_L1 = r_L1 - ri;
            u_L1 = u_L1 - dot(u_L1,u_up_i)*u_up_i;
            u_L1 = u_L1 / max(norm(u_L1),1e-6);

        else
            % ===== 僚机：严格跟踪队形目标点（水平LOS + L1） =====
            r_rel_T = rTi - ri;
            R_to_T = norm(r_rel_T);

            % 队形误差的水平分量
            r_rel_h = r_rel_T - dot(r_rel_T,u_up_i)*u_up_i;
            R_h = max(norm(r_rel_h),1e-6);
            u_to_form = r_rel_h / R_h;

            % 虚拟航段方向：指向队形点
            t_h = u_to_form;

            % 横向误差（与 t_h 正交的水平分量）
            e_h = r_rel_h;
            xt_vec = e_h - dot(e_h, t_h)*t_h;
            xt = norm(xt_vec);

            if xt > 1e-6
                n_xt = xt_vec/xt;
            else
                n_xt = cross(u_up_i,t_h);
                n_xt = n_xt/max(norm(n_xt),1e-6);
            end

            % 僚机 L1 距离：随速度变化但有限制
            L1_dist = max(5e3, 10*Vh);
            L1_dist = min(L1_dist, 40e3);

            % 构造 L1 点：对横向误差强闭环
            r_L1 = ri + L1_dist*t_h - min(K_form_xt*xt, 0.8*L1_dist)*n_xt;

            u_L1 = r_L1 - ri;
            u_L1 = u_L1 - dot(u_L1,u_up_i)*u_up_i;
            u_L1 = u_L1 / max(norm(u_L1),1e-6);

            a_lat_max = a_lat_max_wing;
        end

        % 公共：由 u_L1 得到转向角 eta
        sin_eta = dot(cross(u_vh,u_L1),u_up_i);
        cos_eta = dot(u_vh,u_L1);
        eta = atan2(sin_eta,cos_eta);

        a_lat_cmd = 2*Vh^2/max(L1_dist,1)*sin(eta);
        a_lat_cmd = min(max(a_lat_cmd,-a_lat_max),a_lat_max);

        right_h = cross(u_up_i,u_vh);
        right_h = right_h/max(norm(right_h),1e-6);
        a_lat_vec = a_lat_cmd * right_h;

        %% ============= 高度环（分段PID + 实时重力前馈） =============
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

        [lat_g, lon_g, h_g] = ecef2lla_cgcs2000(ri');
        g_ecef = gravity_cgcs2000_ecef(lat_g, lon_g, h_g);
        g_up = -dot(g_ecef, u_up_i);

        v_up = dot(vi,u_up_i);
        h_err = h_cmd - h;

        int_h(iv) = int_h(iv) + h_err*dt_use;
        int_h(iv) = min(max(int_h(iv),-int_h_lim),int_h_lim);

        a_h_cmd = Kph*h_err + Kih*int_h(iv) - Kdh*v_up + g_up;

        if R_to_T < R_term2
            a_h_lim_now = 1.2*a_h_max;
        else
            a_h_lim_now = a_h_max;
        end
        a_h_cmd = min(max(a_h_cmd,-a_h_lim_now),a_h_lim_now);

        h_floor_200 = 29e3;
        t_protect   = 2;
        if (t_now < t_protect) && (h < h_floor_200)
            a_h_cmd = max(a_h_cmd, 0);
            int_h(iv) = max(int_h(iv), 0);
        end

        a_vert_vec = a_h_cmd * u_up_i;
        a_cmd_ecef = a_lat_vec + a_vert_vec;

        %% ============= 安全屏障：最近邻距离 >= 1km（水平排斥） =============
        if barrier_on
            a_rep = zeros(3,1);
            for jv=1:nVeh
                if jv==iv, continue; end
                rij = ri - r(:,jv);
                rij_h = rij - dot(rij,u_up_i)*u_up_i;
                dij = norm(rij_h);

                if dij < d_safe_min && dij > 1
                    mag = k_rep * (1 - dij/d_safe_min) * a_rep_max;
                    mag = min(max(mag,0), a_rep_max);
                    a_rep = a_rep + mag * (rij_h / dij);
                end
            end
            if norm(a_rep) > a_rep_max
                a_rep = a_rep / norm(a_rep) * a_rep_max;
            end
            a_cmd_ecef = a_cmd_ecef + a_rep;
        end

        %% ============= 解耦：反解 CL_req -> alpha_cmd =============
        a_up_total = dot(a_cmd_ecef, u_up_i);
        a_up_total = max(a_up_total, 0.8*g_up);

        cphi = max(cos(abs(phi(iv))), 0.35);
        L_req  = m(iv) * a_up_total / cphi;
        CL_req = L_req / max(qbar*S_ref(iv), 1.0);

        alpha_cmd = fminbnd(@(a) (aero_coeffs(Ma, a, dT(iv)) - CL_req).^2, alpha_min, alpha_max);
        alpha_sat = min(max(alpha_cmd, alpha_min), alpha_max);

        if abs(alpha_cmd - alpha_sat) > 1e-9
            over = abs(alpha_cmd - alpha_sat)/max(abs(alpha_max-alpha_min),1e-6);
            unload = min(1.0, alpha_unload_gain * (0.2 + 3.0*over));
            a_h_cmd = (1.0 - unload)*a_h_cmd;
            int_h(iv) = int_h(iv) * (1.0 - 0.5*unload);
            alpha_cmd = alpha_sat;
        else
            alpha_cmd = alpha_sat;
        end

        %% ============= a_cmd -> phi/theta/psi（保持原结构） =============
        a_vert_req = dot(a_cmd_ecef,u_up_i);
        a_h_req_vec = a_cmd_ecef - a_vert_req*u_up_i;
        a_lat_req = dot(a_h_req_vec,right_h);

        phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
        phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

        gamma_now = atan2(v_up, Vh);
        theta_cmd = gamma_now + deg2rad(alpha_cmd)*cos(phi_cmd);
        theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);

        chi_v = atan2(dot(u_vh,eE), dot(u_vh,eN));
        psi_cmd = chi_v;

        phi(iv) = phi_cmd;
        theta(iv) = theta_cmd;
        psi(iv) = psi_cmd;

        %% ============= 速度控制（CT 前馈 + PI） =============
        if R_to_T > 900e3
            M_cmd = M_cmd_far;
        elseif R_to_T > 250e3
            M_cmd = M_cmd_mid;
        else
            M_cmd = M_cmd_near;
        end
        V_cmd = max(M_cmd*a_snd, V_floor);

        [CL, CD, ~] = aero_coeffs(Ma, alpha_cmd, dT(iv));
        L = CL*qbar*S_ref(iv);
        D = CD*qbar*S_ref(iv);

        CT_ff = D / max(qbar*S_ref(iv),1);

        v_err = V_cmd - V;
        int_v(iv) = int_v(iv) + v_err*dt_use;
        int_v(iv) = min(max(int_v(iv),-int_v_lim),int_v_lim);

        CT_cmd = CT_ff + Kpv*v_err + Kiv*int_v(iv);
        CT_cmd = min(max(CT_cmd,CT_min),CT_max);

        if (CT_cmd>=CT_max-1e-6 && v_err>0) || (CT_cmd<=CT_min+1e-6 && v_err<0)
            int_v(iv) = int_v(iv) - 0.7*v_err*dt_use;
        end

        CT_grid = zeros(size(dT_grid));
        for ii=1:numel(dT_grid)
            [~,~,CT_grid(ii)] = aero_coeffs(Ma, alpha_cmd, dT_grid(ii));
        end
        [~,ix] = min(abs(CT_grid - CT_cmd));
        dT_target = dT_grid(ix);

        dT(iv) = dT(iv) + (dT_target - dT(iv))*dt_use/tau_dT;
        dT(iv) = min(max(dT(iv),dT_min),dT_max);

        [~,~,CT] = aero_coeffs(Ma, alpha_cmd, dT(iv));
        T_eng = CT * qbar * S_ref(iv);

        %% ============= 力与积分（保持单机3D升力方向含滚转phi） =============
        Vmag = max(norm(vi), 1e-6);
        u_v  = vi / Vmag;

        right = cross(u_up_i, u_v);
        nr = norm(right);
        if nr < 1e-8
            right = cross(u_up_i, eE);
            nr = norm(right);
        end
        right = right / max(nr, 1e-6);

        lift_up = cross(u_v, right);
        lift_up = lift_up / max(norm(lift_up), 1e-6);

        lift_dir = cos(phi(iv)) * lift_up + sin(phi(iv)) * right;
        lift_dir = lift_dir / max(norm(lift_dir), 1e-6);

        F_drag_e   = -D * u_v;
        F_lift_e   =  L * lift_dir;
        F_thrust_e =  T_eng * u_v;

        F_ecef = F_drag_e + F_lift_e + F_thrust_e;

        a_ecef = g_ecef + F_ecef/m(iv) ...
            - 2*cross(omega_ie, vi) ...
            - cross(omega_ie, cross(omega_ie, ri));

        vi = vi + a_ecef*dt_use;
        ri = ri + vi*dt_use;

        v(:,iv) = vi;
        r(:,iv) = ri;

        % 日志
        Rhist(k,:,iv) = ri.';
        Vhist(k,:,iv) = vi.';
        Hhist(k,iv) = h;
        Dist_hist(k,iv) = R_to_T;
        Vmag_hist(k,iv) = V;
        Vcmd_hist(k,iv) = V_cmd;
        dT_hist(k,iv) = dT(iv);
        alpha_hist(k,iv) = alpha_cmd;
    end

    if k > N_max
        break;
    end

    % 时间推进
    t_now = t_now + dt_use;
    t_log(k) = t_now;

    % 记录编队关键距离 & 横向展开
    r1n = r(:,1);
    r2n = r(:,2);
    r3n = r(:,3);
    r4n = r(:,4);

    d12_hist(k) = norm(r2n - r1n);
    d13_hist(k) = norm(r3n - r1n);
    d24_hist(k) = norm(r4n - r2n);
    d34_hist(k) = norm(r4n - r3n);

    y2 = dot(r2n - r1n, ey);
    y3 = dot(r3n - r1n, ey);
    width_hist(k) = abs(y3 - y2);

    k = k + 1;
end

%% ================= 截断有效数据 =================
valid = isfinite(t_log);
tt = t_log(valid);

Hk = Hhist(valid,:);
Dk = Dist_hist(valid,:);

idx = find(valid);
d12 = d12_hist(idx);
d13 = d13_hist(idx);
d24 = d24_hist(idx);
d34 = d34_hist(idx);
width = width_hist(idx);

fprintf('Simulation stop reason: %s, t=%.2f s, leader min range=%.1f m\n', ...
    stop_reason, tt(end), min(Dk(:,1)));

fprintf('Formation check (min over time):\n');
fprintf('  min d12=%.1f m, min d13=%.1f m, min d24=%.1f m, min d34=%.1f m\n', ...
    min(d12), min(d13), min(d24), min(d34));
fprintf('  max width=%.1f m (limit 2000m)\n', max(width));

%% ================= 绘图：领机距离/高度、编队距离与宽度 =================
figure;
plot(tt, Dk(:,1)/1e3,'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('Leader range to target (km)');
title('领机到目标距离');

figure;
plot(tt, Hk(:,1)/1e3,'LineWidth',1.3); grid on; yline(30,'r--');
xlabel('Time (s)'); ylabel('Leader altitude (km)');
title('领机高度');

figure;
plot(tt, d12/1e3,'LineWidth',1.2); hold on; grid on;
plot(tt, d13/1e3,'LineWidth',1.2);
plot(tt, d24/1e3,'LineWidth',1.2);
plot(tt, d34/1e3,'LineWidth',1.2);
yline(1.0,'k--','d_{safe}=1km');
xlabel('Time (s)'); ylabel('Distance (km)');
title('编队关键邻接距离');
legend('d(1-2)','d(1-3)','d(2-4)','d(3-4)','Location','best');

figure;
plot(tt, width/1e3,'LineWidth',1.3); grid on;
yline(2.0,'k--','width<=2km');
xlabel('Time (s)'); ylabel('Wing-to-wing width (km)');
title('横向展开约束检查');

%% ================= 所有飞行器：轨迹经纬高图（LLA, Lon-Lat-Alt） =================
Rk = Rhist(valid,:,:);
Ntraj = size(Rk,1);

lat_traj = nan(Ntraj, nVeh);
lon_traj = nan(Ntraj, nVeh);
h_traj   = nan(Ntraj, nVeh);

for iv = 1:nVeh
    for ii = 1:Ntraj
        [lat_traj(ii,iv), lon_traj(ii,iv), h_traj(ii,iv)] = ...
            ecef2lla_cgcs2000( squeeze(Rk(ii,:,iv)) );
    end
end

figure; hold on; grid on;
clr = lines(nVeh);

for iv = 1:nVeh
    plot3(lon_traj(:,iv), lat_traj(:,iv), h_traj(:,iv)/1000, ...
        'Color', clr(iv,:), 'LineWidth', 1.5);
end

plot3(lon0, lat0, h0/1000, 'ko', 'MarkerFaceColor','g', 'MarkerSize',7);
plot3(lonT, latT, hT/1000, 'ko', 'MarkerFaceColor','r', 'MarkerSize',7);

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
zlabel('Altitude (km)');
title('All Vehicles Trajectories in LLA (Lon-Lat-Alt)');

leg = cell(nVeh+2,1);
for iv = 1:nVeh
    leg{iv} = sprintf('Veh %d', iv);
end
leg{nVeh+1} = 'Launch';
leg{nVeh+2} = 'Target';
legend(leg, 'Location','best');

view(3);