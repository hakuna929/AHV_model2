% test3.m
% 4-ship formation (diamond, vertical split) + leader to target
% ------------------------------------------------------------
% 设定：
% - 3DOF 点质量（ECEF），姿态(phi/theta/psi)视为直接可控
% - 领机(1号)飞向目标 rT（保持原"末段增强"逻辑）
% - 僚机(2~4号)跟踪领机航迹坐标系下的菱形偏移（含高度错层）
% - 攻角严格限制 alpha ∈ [-2°, +1°]
% - 推力模型移除调试系数5：T = CT*qbar*S_ref
%
% 依赖函数：
% lla2ecef_cgcs2000, ecef2lla_cgcs2000, atmos_simple, aero_coeffs
%
% 注意：
% - lla2ecef/ecef2lla 的 lat/lon 单位为度，因此三角函数用 sind/cosd
% - alpha/phi/theta 等弧度

clear; clc;

%% ===================== 仿真时间 =====================
dt_base = 0.01;
T_end   = 3000;
N_max   = ceil(T_end/dt_base) * 5;

%% ===================== 常量 =====================
mu  = 3.986004418e14;
Re  = 6378137.0;
we  = 7.2921150e-5;
omega_ie = [0;0;we];

use_simple_gravity = true;
use_rotation_terms = false;
g0 = 9.80665;

%% ===================== 飞行器参数 =====================
m     = 671.33;
S_ref = 0.2986;

%% ===================== 初始/目标 =====================
h0 = 30e3; lat0 = 19.2; lon0 = 110.5;   % deg
r0 = lla2ecef_cgcs2000(lat0, lon0, h0); r0 = r0(:);

latT = 13.35; lonT = 144.55; hT = 30e3;
rT = lla2ecef_cgcs2000(latT, lonT, hT); rT = rT(:);

[a0,~] = atmos_simple(h0);
V0 = 6.5 * a0;

% 航点（领机可扩展）
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
u_up0 = r0/norm(r0);
seg0 = r_wp_next - r_wp_prev;
seg0_h = seg0 - dot(seg0,u_up0)*u_up0;
if norm(seg0_h) < 1e-6, seg0_h = seg0; end
v0 = V0 * seg0_h / norm(seg0_h);

%% ===================== 编队参数（菱形 + 高度错层） =====================
Nveh = 4;
d  = 3000;   % 平面间距 3 km
du = 500;    % 垂向错层 500 m

% Delta = [前, 右, 上]（航迹坐标）
Delta = [   0,    0,   0;     % 1: lead
           -d,   -d, +du;     % 2: left-back (up)
           -d,   +d, -du;     % 3: right-back (down)
         -2*d,    0,   0];    % 4: back

% 初始把编队按 v0 定义的航迹坐标系散开
ef0 = seg0_h / max(norm(seg0_h),1e-6);
eu0 = u_up0;
er0 = cross(eu0,ef0); er0 = er0 / max(norm(er0),1e-6);

r = zeros(3,Nveh);
v = zeros(3,Nveh);
for i=1:Nveh
    r(:,i) = r0 + Delta(i,1)*ef0 + Delta(i,2)*er0 + Delta(i,3)*eu0;
    v(:,i) = v0;
end

%% ===================== 控制参数 =====================
% ---------- 领机横向：L1 ----------
L1_base  = 100e3;
L1_gainV = 25;
a_lat_max_cruise = 60;
a_lat_max_term   = 140;

% ---------- 僚机相对控制（BB：不用L1航线，直接相对PD） ----------
% 在航迹坐标系中，对 (右向误差、上向误差) 给加速度指令
Kp_lat_rel = 0.010;   % (m/s^2)/m  侧向"弹簧"
Kd_lat_rel = 0.080;   % (m/s^2)/(m/s) 侧向"阻尼"

Kp_up_rel  = 0.010;   % 垂向
Kd_up_rel  = 0.080;

a_lat_rel_max = 40;   % 僚机相对横向加速度限幅（建议小于领机末段）
a_up_rel_max  = 18;   % 垂向限幅

% ---------- 高度外环（领机用，僚机用ref高度） ----------
Kph = 0.010; Kih = 0.00002; Kdh = 0.008;
int_h = zeros(Nveh,1);
int_h_lim = 8e4;
a_h_max_lead = 18;

% alpha/phi 限幅（弧度）
phi_lim   = 65*pi/180;
theta_lim = 35*pi/180;
alpha_min = -2*pi/180;
alpha_max =  1*pi/180;

alpha_unload_gain = 0.85;

% ---------- 速度控制（CT前馈+PI） ----------
Kpv = 0.0018; Kiv = 0.00025;
int_v = zeros(Nveh,1); 
int_v_lim = 8e3;

CT_min = 0.05; CT_max = 0.1;
dT_min = 0.10; dT_max = 1.00;
tau_dT = 0.3; 
dT = 0.5*ones(Nveh,1);

% 僚机速度"前向间距补偿"
Kf = 0.02;          % 1/s
dV_max = 100;       % m/s 限幅

% 马赫指令（领机/僚机都用同一基准，再加dV）
M_cmd_far  = 6.5;
M_cmd_mid  = 6.5;
M_cmd_near = 6.5;
V_floor    = 1100;

% ---------- 终端捕获参数（领机判断） ----------
R_hit            = 1000;
R_term1          = 120e3;
R_term2          = 40e3;
R_pass_check     = 20e3;
R_abort_diverge  = 50e3;

%% ===================== 记录数组（领机为主 + 队形误差） =====================
t_log = nan(N_max,1);

Rhist1 = nan(N_max,3);
Vhist1 = nan(N_max,3);
Hhist1 = nan(N_max,1);
Dist_hist = nan(N_max,1);

% 编队误差（僚机相对ref的前/右/上误差）
ef_err = nan(N_max, Nveh);
er_err = nan(N_max, Nveh);
eu_err = nan(N_max, Nveh);

alpha_hist = nan(N_max, Nveh);
phi_hist   = nan(N_max, Nveh);
theta_hist = nan(N_max, Nveh);
Vmag_hist  = nan(N_max, Nveh);
Vcmd_hist  = nan(N_max, Nveh);
dT_hist    = nan(N_max, Nveh);

stop_reason = "completed";
k = 1;
t_now = 0;

prev_R = inf;
min_R = inf;

%% ===================== 主循环 =====================
while k <= N_max && t_now <= T_end

    % ===== 先用领机(1号)算当前地理量/环境（用于编队坐标系）=====
    rL = r(:,1);
    vL = v(:,1);

    [lat_deg_L, lon_deg_L, hL] = ecef2lla_cgcs2000(rL');
    if hL <= 20e3
        stop_reason = "terrain impact";
        fprintf('Leader impacted terrain at t=%.1f s\n', t_now);
        break;
    end

    % 目标相对（领机）
    r_rel_T = rT - rL;
    R_to_T = norm(r_rel_T);
    u_los_T = r_rel_T / max(R_to_T,1);
    Vc_toward = dot(vL, u_los_T);

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

    % --- 自适应步长（按领机终端距离）---
    if R_to_T < R_term2
        dt_use = 0.001;
    elseif R_to_T < R_term1
        dt_use = 0.01;
    else
        dt_use = dt_base;
    end

    % --- 环境（对每架机分别算也可以；这里每架机单独算）---

    % ====== 领机航迹坐标系（前/右/上）======
    eu = rL / norm(rL);
    vL_h = vL - dot(vL,eu)*eu;
    ef = vL_h / max(norm(vL_h),1e-6);
    er = cross(eu,ef); er = er / max(norm(er),1e-6);

    % ====== 生成 4 架参考位置 r_ref（ECEF）======
    r_ref = zeros(3,Nveh);
    for i=1:Nveh
        r_ref(:,i) = rL + Delta(i,1)*ef + Delta(i,2)*er + Delta(i,3)*eu;
    end

    % ====== 航点切换（领机）======
    if norm(r_wp_next - rL) < wp_switch_dist && wp_idx < nWP
        wp_idx = wp_idx + 1;
        r_wp_prev = r_wp_next;
        r_wp_next = wps_ecef(wp_idx,:)';
    end

    % ====== 对每架机计算控制并积分 ======
    for i = 1:Nveh

        ri = r(:,i);
        vi = v(:,i);

        % 当前高度/经纬（用于构造 eE/eN）
        [lat_deg, lon_deg, h] = ecef2lla_cgcs2000(ri');

        % 大气
        V = max(norm(vi),1e-3);
        [a_snd, rho] = atmos_simple(max(h,0));
        Ma = max(min(V/max(a_snd,1e-3),8.0),0.0);
        qbar = 0.5*rho*V^2;

        % 局部基（每架机自己的 up）
        u_up = ri / norm(ri);

        % EN（用于psi显示/兼容原逻辑）
        eE = [-sind(lon_deg); cosd(lon_deg); 0];
        eN = [-sind(lat_deg)*cosd(lon_deg); -sind(lat_deg)*sind(lon_deg); cosd(lat_deg)];

        v_h = vi - dot(vi,u_up)*u_up;
        Vh = max(norm(v_h),1e-6);
        u_vh = v_h / Vh;

        right_h = cross(u_up,u_vh); right_h = right_h/max(norm(right_h),1e-6);

        % ---------------- 横向/纵向加速度指令 ----------------
        if i == 1
            % ========== 领机：L1 航线跟踪 ==========
            seg = r_wp_next - r_wp_prev;
            t_hat = seg / max(norm(seg),1e-6);
            t_h = t_hat - dot(t_hat,u_up)*u_up;
            t_h = t_h / max(norm(t_h),1e-6);

            r_from_start = ri - r_wp_prev;
            r_xt = r_from_start - dot(r_from_start,t_h)*t_h;
            xt = norm(r_xt);

            if xt > 1e-6
                n_xt = r_xt/xt;
            else
                n_xt = cross(u_up,t_h); n_xt = n_xt/max(norm(n_xt),1e-6);
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

            r_L1 = ri + L1_dist*t_h - min(xt,0.5*L1_dist)*n_xt;
            u_L1 = r_L1 - ri;
            u_L1 = u_L1 - dot(u_L1,u_up)*u_up;
            u_L1 = u_L1 / max(norm(u_L1),1e-6);

            sin_eta = dot(cross(u_vh,u_L1),u_up);
            cos_eta = dot(u_vh,u_L1);
            eta = atan2(sin_eta,cos_eta);

            a_lat_cmd = 2*Vh^2/max(L1_dist,1)*sin(eta);
            a_lat_cmd = min(max(a_lat_cmd,-a_lat_max),a_lat_max);

            a_lat_vec = a_lat_cmd * right_h;

            % 领机高度控制：保持目标高度 hT（也可改30km常数）
            h_cmd = hT;
            v_up = dot(vi,u_up);
            h_err = h_cmd - h;

            int_h(i) = int_h(i) + h_err*dt_use;
            int_h(i) = min(max(int_h(i),-int_h_lim),int_h_lim);

            a_h_cmd = Kph*h_err + Kih*int_h(i) - Kdh*v_up;
            a_h_cmd = min(max(a_h_cmd,-a_h_max_lead),a_h_max_lead);

            a_vert_vec = a_h_cmd * u_up;

        else
            % ========== 僚机：相对PD（BB：直接相对误差->加速度）==========
            % 误差在"领机航迹坐标系"(ef/er/eu)上分解
            e = r_ref(:,i) - ri;

            % 相对速度误差（也投影到同一基）
            ev = v(:,1) - vi;

            e_f = dot(e, ef);
            e_r = dot(e, er);
            e_u = dot(e, eu);

            ev_f = dot(ev, ef);
            ev_r = dot(ev, er);
            ev_u = dot(ev, eu);

            % 侧向/垂向加速度（PD）
            a_r_cmd = Kp_lat_rel*e_r + Kd_lat_rel*ev_r;
            a_u_cmd = Kp_up_rel *e_u + Kd_up_rel *ev_u;

            a_r_cmd = min(max(a_r_cmd, -a_lat_rel_max), a_lat_rel_max);
            a_u_cmd = min(max(a_u_cmd, -a_up_rel_max ), a_up_rel_max );

            % 转回ECEF加速度向量（用领机基向量 er/eu）
            a_lat_vec = a_r_cmd * er;

            a_vert_vec = a_u_cmd * eu;

            % 僚机高度"参考"是 r_ref 的高度（BB）
            % 这里不再用高度环，直接用 a_u_cmd 做垂向控制（最简）
            % 如果你想更平滑，可把 a_u_cmd 换成对(h_ref-h)的PID再限幅。

        end

        % 合成加速度指令
        a_cmd_ecef = a_lat_vec + a_vert_vec;
        a_cmd_norm = norm(a_cmd_ecef);

        % ---------------- a_cmd -> alpha/phi（含alpha饱和卸载） ----------------
        % 先算 CL_req
        CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);

        % 线性反推 alpha（近似）
        adeg_est = 0.0;
        CL_alpha_rad = (0.07235 - 0.003368*Ma) * (180/pi);
        CL_alpha_rad = max(CL_alpha_rad,0.1);
        CL0_est = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.001185*adeg_est^2;

        alpha_cmd = (CL_req - CL0_est)/CL_alpha_rad;
        alpha_sat = min(max(alpha_cmd, alpha_min), alpha_max);

        if abs(alpha_cmd - alpha_sat) > 1e-9
            % 卸载"垂向"分量（尽量保队形/横向）
            over = abs(alpha_cmd - alpha_sat)/max(abs(alpha_max-alpha_min),1e-6);
            unload = min(1.0, alpha_unload_gain * (0.2 + 3.0*over));

            % 只缩小垂向加速度（a_vert_vec），让alpha回到范围
            a_vert_vec = (1.0 - unload)*a_vert_vec;

            if i==1
                % 抗积分风up
                int_h(i) = int_h(i) * (1.0 - 0.5*unload);
            end

            a_cmd_ecef = a_lat_vec + a_vert_vec;
            a_cmd_norm = norm(a_cmd_ecef);

            CL_req = (m*a_cmd_norm)/max(qbar*S_ref,1);
            alpha_cmd = (CL_req - CL0_est)/CL_alpha_rad;
            alpha_cmd = min(max(alpha_cmd, alpha_min), alpha_max);
        else
            alpha_cmd = alpha_sat;
        end

        % 滚转由"加速度分量"分配（近似）
        a_vert_req = dot(a_cmd_ecef,u_up);
        a_h_req_vec = a_cmd_ecef - a_vert_req*u_up;
        a_lat_req = dot(a_h_req_vec,right_h);

        phi_cmd = atan2(a_lat_req, max(abs(a_vert_req),1.0));
        phi_cmd = min(max(phi_cmd,-phi_lim),phi_lim);

        v_up = dot(vi,u_up);
        gamma_now = atan2(v_up, Vh);
        theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
        theta_cmd = min(max(theta_cmd,-theta_lim),theta_lim);

        chi_v = atan2(dot(u_vh,eE), dot(u_vh,eN));
        psi_cmd = chi_v;

        % ---------------- 速度指令 ----------------
        if R_to_T > 900e3
            M_cmd = M_cmd_far;
        elseif R_to_T > 250e3
            M_cmd = M_cmd_mid;
        else
            M_cmd = M_cmd_near;
        end
        V_base = max(M_cmd*a_snd, V_floor);

        if i==1
            V_cmd = V_base;
        else
            % 前向间距补偿（用领机航迹基向量 ef）
            e = r_ref(:,i) - ri;
            e_f = dot(e, ef);
            dV = Kf * e_f;
            dV = min(max(dV, -dV_max), dV_max);
            V_cmd = V_base + dV;
        end

        % ---------------- 速度环：CT前馈+PI -> dT ----------------
        [CL, CD, CT_now] = aero_coeffs(Ma, alpha_cmd, dT(i));

        L = CL*qbar*S_ref;
        D = CD*qbar*S_ref;

        CT_ff = D / max(qbar*S_ref,1);

        v_err = V_cmd - V;
        int_v(i) = int_v(i) + v_err*dt_use;
        int_v(i) = min(max(int_v(i),-int_v_lim),int_v_lim);

        CT_cmd = CT_ff + Kpv*v_err + Kiv*int_v(i);
        CT_cmd = min(max(CT_cmd,CT_min),CT_max);

        if (CT_cmd>=CT_max-1e-6 && v_err>0) || (CT_cmd<=CT_min+1e-6 && v_err<0)
            int_v(i) = int_v(i) - 0.7*v_err*dt_use;
        end

        dT_grid = linspace(dT_min,dT_max,41);
        CT_grid = zeros(size(dT_grid));
        for ii=1:numel(dT_grid)
            [~,~,CT_grid(ii)] = aero_coeffs(Ma, alpha_cmd, dT_grid(ii));
        end
        [~,ix] = min(abs(CT_grid - CT_cmd));
        dT_target = dT_grid(ix);

        dT(i) = dT(i) + (dT_target - dT(i))*dt_use/tau_dT;
        dT(i) = min(max(dT(i),dT_min),dT_max);

        [~,~,CT] = aero_coeffs(Ma, alpha_cmd, dT(i));

        % 推力（去掉系数5）
        T_eng = CT * qbar * S_ref;

        % ---------------- 力与积分 ----------------
        if norm(v_h) > 1e-6
            fwd = v_h / norm(v_h);
        else
            fwd = vi / max(norm(vi),1e-6);
        end

        F_drag_e   = -D * (vi/max(norm(vi),1e-6));
        F_lift_e   =  L * u_up;
        F_thrust_e =  T_eng * fwd;
        F_ecef = F_drag_e + F_lift_e + F_thrust_e;

        if use_simple_gravity
            g_ecef = -mu/norm(ri)^3 * ri;
        end

        if use_rotation_terms
            a_ecef = g_ecef + F_ecef/m - 2*cross(omega_ie,vi) - cross(omega_ie,cross(omega_ie,ri));
        else
            a_ecef = g_ecef + F_ecef/m;
        end

        % 积分更新
        v(:,i) = vi + a_ecef*dt_use;
        r(:,i) = ri + v(:,i)*dt_use;

        % 记录一些量
        alpha_hist(k,i) = alpha_cmd;
        phi_hist(k,i)   = phi_cmd;
        theta_hist(k,i) = theta_cmd;
        Vmag_hist(k,i)  = V;
        Vcmd_hist(k,i)  = V_cmd;
        dT_hist(k,i)    = dT(i);

        if any(~isfinite([r(:,i);v(:,i);Ma;qbar;CL;CD;CT;dT(i)]))
            stop_reason = "NaN/Inf";
            fprintf('NaN/Inf at t=%.2f s (veh %d)\n', t_now, i);
            k = k + 1;
            break;
        end

    end % for each vehicle

    if stop_reason == "NaN/Inf"
        break;
    end

    % 时间推进
    t_now = t_now + dt_use;

    % -------- 记录（以领机为主）--------
    t_log(k) = t_now;
    Rhist1(k,:) = r(:,1).';
    Vhist1(k,:) = v(:,1).';
    Hhist1(k) = hL;
    Dist_hist(k) = R_to_T;

    % 编队误差（投影到领机航迹基）
    for i=1:Nveh
        e = r_ref(:,i) - r(:,i);
        ef_err(k,i) = dot(e, ef);
        er_err(k,i) = dot(e, er);
        eu_err(k,i) = dot(e, eu);
    end

    if mod(round(t_now,1),1.0)==0
        fprintf("t=%.0f  R_to_T=%.0fkm  hL=%.0f  leadV=%.0f  e2(f,r,u)=(%.0f,%.0f,%.0f)m\n", ...
            t_now, R_to_T/1000, hL, norm(v(:,1)), ef_err(k,2), er_err(k,2), eu_err(k,2));
    end

    k = k + 1;

end % while

%% ===================== 截断有效数据 =====================
valid = isfinite(t_log) & isfinite(Dist_hist);
tt = t_log(valid);

R1 = Rhist1(valid,:);
V1 = Vhist1(valid,:);
H1 = Hhist1(valid);
Dk = Dist_hist(valid);

efE = ef_err(valid,:);
erE = er_err(valid,:);
euE = eu_err(valid,:);

alphak = alpha_hist(valid,:);
Vmk    = Vmag_hist(valid,:);
Vcmdk  = Vcmd_hist(valid,:);

fprintf('Simulation stop reason: %s, t=%.2f s, min range=%.1f m\n', stop_reason, tt(end), min(Dk));

%% ===================== 绘图 =====================
figure; plot(tt, Dk/1000, 'm', 'LineWidth',1.7); grid on;
xlabel('Time (s)'); ylabel('Leader distance to Target (km)'); title('Leader Distance');

figure; plot(tt, H1/1000, 'LineWidth',1.4); grid on;
xlabel('Time (s)'); ylabel('Leader altitude (km)'); title('Leader Altitude');

figure; 
plot(tt, efE(:,2), 'LineWidth',1.2); hold on; grid on;
plot(tt, erE(:,2), 'LineWidth',1.2);
plot(tt, euE(:,2), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Error (m)'); title('Wingman 2 formation error (proj on leader frame)');
legend('e_f','e_r','e_u');

figure;
plot(tt, efE(:,3), 'LineWidth',1.2); hold on; grid on;
plot(tt, erE(:,3), 'LineWidth',1.2);
plot(tt, euE(:,3), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Error (m)'); title('Wingman 3 formation error');
legend('e_f','e_r','e_u');

figure;
plot(tt, efE(:,4), 'LineWidth',1.2); hold on; grid on;
plot(tt, erE(:,4), 'LineWidth',1.2);
plot(tt, euE(:,4), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Error (m)'); title('Wingman 4 formation error');
legend('e_f','e_r','e_u');

figure;
plot(tt, alphak*180/pi, 'LineWidth',1.1); grid on;
yline(-2,'r--'); yline(1,'r--');
xlabel('Time (s)'); ylabel('\alpha (deg)'); title('Angle of Attack (all vehicles)');
legend('1','2','3','4','Location','best');

figure;
plot(tt, Vmk(:,1), 'k', 'LineWidth',1.2); hold on; grid on;
plot(tt, Vcmdk(:,1), 'r--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (m/s)'); title('Leader speed tracking');
legend('|V|','V_{cmd}','Location','best');

figure;
plot3(R1(:,1)/1e3, R1(:,2)/1e3, R1(:,3)/1e3, 'b','LineWidth',1.3); hold on; grid on; axis equal;
[xe,ye,ze] = sphere(60);
surf(Re*xe/1e3, Re*ye/1e3, Re*ze/1e3, 'FaceAlpha',0.08,'EdgeColor','none');
plot3(rT(1)/1e3, rT(2)/1e3, rT(3)/1e3, 'ro','MarkerFaceColor','r');
xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
title('Leader Trajectory in ECEF');
legend('Leader','Earth','Target','Location','best');