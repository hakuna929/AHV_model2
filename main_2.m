clear;
clc;

%% 时间
dt = 0.1;
T  = 3000;
t  = 0:dt:T;
N  = length(t);

%% 地球常量
mu  = 3.986004418e14;
J2  = 1.08263e-3;
Re  = 6378137.0;            % m

% 地球自转参数
we  = 7.2921150e-5;
omega_ie = [0;0;we];

%% 飞行器参数 
m     = 671.33;
S_ref = 0.2986;

%% 目标点设置（经纬高）
h_target = 30e3;
lat_T = 13.35;
lon_T = 144.55;
rT = lla2ecef_cgcs2000(lat_T, lon_T, h_target); rT = rT(:);
vT = [0;0;0];

%% 初始条件
h_init = 30e3;
lat0 = 19.2;
lon0 = 110.5;
r = lla2ecef_cgcs2000(lat0, lon0, h_init);
r0_ecef = r;

Ma_r = 6.5;
[a0,~] = atmos_simple(h_init);
V0 = Ma_r * a0;
V_cmd = V0;

r0_hat = r0_ecef / norm(r0_ecef);
v_dir = rT - (rT' * r0_hat) * r0_hat;
v_dir = v_dir / norm(v_dir);
v = V0 * v_dir;

% 姿态初始化
eU = [cosd(lat0)*cosd(lon0); cosd(lat0)*sind(lon0); sind(lat0)];
eN = [-sind(lat0)*cosd(lon0); -sind(lat0)*sind(lon0); cosd(lat0)];
r_rel0 = rT - r0_ecef;
x_h = r_rel0 - (r_rel0' * r0_hat) * r0_hat;
if norm(x_h) < 1e-6
    xb = eN;
else
    xb = x_h / norm(x_h);
end
yb = cross(eU, xb); yb = yb / max(norm(yb), 1e-6);
zb = cross(xb, yb); zb = zb / max(norm(zb), 1e-6);
C_b2e_0 = [xb, yb, zb];

theta = asin(-C_b2e_0(3,1));
psi   = atan2(C_b2e_0(2,1), C_b2e_0(1,1));
phi   = atan2(C_b2e_0(3,2), C_b2e_0(3,3));

% 舵面默认
de1=0.1; de2=0.1; dr=0.1; dT=0.5;

%% 控制参数
Npn_far = 3.0;
Kv = 0.0020;                % 速度环稍增强
a_cmd_max = 120;            % 关键：提高机动上限，先保证不U型发散

% 航向控制
k_psi = 2.0;
psi_rate_max = 10*pi/180;   % rad/s

% 末段发散保护参数
diverge_count = 0;
diverge_count_th = 30;      % 连续3秒 (dt=0.1)
boost_pp_mode = false;

%% 记录
Rhist = zeros(N,3);
Vhist = zeros(N,3);
Hhist = zeros(N,1);
Dist_hist = zeros(N,1);
Vc_hist = zeros(N,1);
Npn_use_hist = zeros(N,1);
chi_v_hist = zeros(N,1);
chi_los_hist = zeros(N,1);
psi_hist = zeros(N,1);
phi_cmd_hist = zeros(N,1);
theta_cmd_hist = zeros(N,1);
alpha_hist = zeros(N,1);
alpha_cmd_hist = zeros(N,1);
CL_req_hist = zeros(N,1);
CL_hist = zeros(N,1);
a_up_aero_hist = zeros(N,1);
Vcmd_hist = zeros(N,1);
mode_hist = zeros(N,1);     % 0:PN主导 1:PP主导

inf_time = N;

for k = 1:N
    %% 地理量
    [lat, lon, h] = ecef2lla_cgcs2000(r');
    Hhist(k)=h;

    if h <= 10e3
        fprintf('Vehicle crashed into the ground at t=%.2f s\n', t(k));
        inf_time = k;
        break;
    end

    %% 空速
    v_air = v;
    Vair = max(norm(v_air),1e-3);

    %% 姿态矩阵
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

    %% alpha beta
    v_air_b = Cb_e * v_air;
    u=v_air_b(1); vv=v_air_b(2); w=v_air_b(3);
    alpha = atan2(w,u);
    beta  = asin(max(min(vv/Vair,1),-1));
    alpha = max(min(alpha, 8*pi/180), -5*pi/180);
    beta  = max(min(beta,  5*pi/180), -5*pi/180);

    %% 大气
    [a,rho] = atmos_simple(max(h,0));
    Ma = max(min(Vair/max(a,1e-3), 7.0), 0.0);
    qbar = 0.5*rho*Vair^2;

    %% 分段减速（关键：减小转弯半径）
    if norm(rT-r) > 800e3
        V_cmd = 6.5 * a;
    elseif norm(rT-r) > 300e3
        V_cmd = 4.0 * a;
    elseif norm(rT-r) > 120e3
        V_cmd = 2.5 * a;
    else
        V_cmd = 1.8 * a;
    end

    %% 重力
    x=r(1); y=r(2); z=r(3); rr=norm(r);
    g0 = -mu/rr^3*[x;y;z];
    gJ2 = [3*mu*J2*Re^2/(2*rr^5)*x*(1-5*(z/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*y*(1-5*(z/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*z*(3-5*(z/rr)^2)];
    g_ecef = g0 + gJ2;

    %% LOS / PN
    r_rel = rT - r;
    R_dist = max(norm(r_rel),1);
    u_los = r_rel / R_dist;
    v_rel = vT - v;
    Vc = -dot(r_rel,v_rel)/R_dist;
    omega_los = cross(r_rel,v_rel)/(R_dist^2 + 1e-6);

    % 发散检测：连续Vc<0则认为飞过最近点
    if Vc < 0
        diverge_count = diverge_count + 1;
    else
        diverge_count = max(diverge_count-1,0);
    end
    if diverge_count >= diverge_count_th
        boost_pp_mode = true;
    end
    if R_dist < 80e3
        boost_pp_mode = true; % 末段直接PP主导
    end

    % PN增益调度
    if R_dist < 20e3
        Npn_use = 1.2;
    elseif R_dist < 200e3
        Npn_use = 1.8;
    elseif R_dist < 800e3
        Npn_use = 2.5;
    else
        Npn_use = Npn_far;
    end

    % 关键：叉乘顺序
    a_pn_ecef = Npn_use * max(Vc,0) * cross(u_los, omega_los);

    if norm(a_pn_ecef) > a_cmd_max
        a_pn_ecef = a_pn_ecef/norm(a_pn_ecef) * a_cmd_max;
    end

    % 补偿项
    a_req_ecef = a_pn_ecef - g_ecef + 2*cross(omega_ie, v) + cross(omega_ie, cross(omega_ie, r));

    %% 当地基底
    u_up = r/norm(r);
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

    % LOS水平航向
    los_h = u_los - dot(u_los,u_up)*u_up;
    if norm(los_h) > 1e-8
        los_h_hat = los_h / norm(los_h);
        chi_los = atan2(dot(los_h_hat,eE_now), dot(los_h_hat,eN_now));
    else
        chi_los = chi_v;
    end

    %% 垂直/水平分解 + 能力分配
    a_vert_geo = dot(a_req_ecef, u_up);
    a_h_geo = a_req_ecef - a_vert_geo*u_up;
    a_h = norm(a_h_geo);

    % 定高外环
    dh = h_init - h;
    vh = dot(v_air, u_up);
    Kph = 0.0012;
    Kdh = 0.08;
    a_hold = Kph*dh - Kdh*vh;
    a_vert_cmd = a_vert_geo + a_hold;

    % 能力上限
    CL_max = 2.5;
    a_L_cap = (CL_max*qbar*S_ref)/m;

    % 末段/发散时优先水平机动
    if boost_pp_mode
        vert_share = 0.40;
    elseif R_dist < 300e3
        vert_share = 0.50;
    else
        vert_share = 0.65;
    end
    a_vert_cmd = max(min(a_vert_cmd, vert_share*a_L_cap), -vert_share*a_L_cap);

    a_h_cap = sqrt(max(a_L_cap^2 - a_vert_cmd^2, 0));
    if a_h > a_h_cap && a_h > 1e-6
        a_h_geo = a_h_geo * (a_h_cap/a_h);
    end

    %% 航向控制：PN主导 <-> PP主导切换
    if boost_pp_mode
        % PP主导：直接追LOS，尽快避免继续远离
        psi_ref = chi_los;
        mode_hist(k) = 1;
    else
        % PN主导：速度航向+LOS误差引导
        psi_ref = wrapToPi(chi_v + 0.9*wrapToPi(chi_los - chi_v));
        mode_hist(k) = 0;
    end

    e_psi = wrapToPi(psi_ref - psi);
    psi_rate_cmd = max(min(k_psi*e_psi, psi_rate_max), -psi_rate_max);
    psi_cmd = wrapToPi(psi + psi_rate_cmd*dt);

    % 滚转由水平加速度确定
    right_dir = cross(u_up, v_h_hat);
    right_dir = right_dir / max(norm(right_dir),1e-6);
    a_lat = dot(a_h_geo, right_dir);

    phi_cmd = atan2(a_lat, max(abs(a_vert_cmd), 2.0));
    phi_lim = 70*pi/180;
    phi_cmd = max(min(phi_cmd, phi_lim), -phi_lim);

    %% 迎角与俯仰
    a_h_now = norm(a_h_geo);
    a_L_req = sqrt(a_vert_cmd^2 + a_h_now^2);
    CL_req = (m*a_L_req)/max(qbar*S_ref,1);

    de1_deg = de1*180/pi; de2_deg = de2*180/pi;
    CL0 = 0.1498 - 0.02751*Ma + 0.002343*Ma^2 + 0.006515*(de1_deg + de2_deg);
    CL_alpha_rad = (0.07235 - 0.003368*Ma)*(180/pi);
    CL_alpha_rad = max(CL_alpha_rad, 0.1);

    alpha_cmd = (CL_req - CL0)/CL_alpha_rad;
    alpha_cmd = max(min(alpha_cmd, 18*pi/180), -3*pi/180);

    v_vert = dot(v_air,u_up);
    v_horz = norm(v_air - v_vert*u_up);
    gamma_now = atan2(v_vert, v_horz);

    theta_cmd = gamma_now + alpha_cmd*cos(phi_cmd);
    theta_cmd = max(min(theta_cmd, 35*pi/180), -35*pi/180);

    %% 理想姿态跟踪
    phi = phi_cmd;
    theta = theta_cmd;
    psi = psi_cmd;

    %% 速度环
    err_v = V_cmd - Vair;
    dT = dT + Kv*err_v*dt;
    dT = max(min(dT,1.0),0.05);

    %% 气动/推力
    alpha_aero = alpha_cmd;
    beta_aero = beta;

    [CL,CD,CY,~,~,~] = aero_coeffs(Ma,alpha_aero,beta_aero,de1,de2,dr,dT);
    L = CL*qbar*S_ref; D=CD*qbar*S_ref; Y=CY*qbar*S_ref;
    Fw = [-D;Y;-L];

    ca=cos(alpha_aero); sa=sin(alpha_aero); cb=cos(beta_aero); sb=sin(beta_aero);
    Cb_w = [ ca*cb, -ca*sb, -sa;
             sb,     cb,     0;
             sa*cb, -sa*sb,  ca ];
    Fa_b = Cb_w*Fw;

    CT = thrust_coeffs(Ma,alpha_aero,beta_aero,dT);
    Ft_b = [CT*qbar*S_ref;0;0];

    F_ecef = Ce_b*(Fa_b + Ft_b);
    a_aero_ecef = (Ce_b*Fa_b)/m;
    a_up_aero = dot(a_aero_ecef, u_up);

    %% 平动
    a_ecef = g_ecef + F_ecef/m - 2*cross(omega_ie,v) - cross(omega_ie,cross(omega_ie,r));
    v = v + a_ecef*dt;
    r = r + v*dt;

    %% 记录
    Rhist(k,:) = r.';
    Vhist(k,:) = v.';
    Dist_hist(k) = norm(rT-r);
    Hhist(k) = h;
    Vc_hist(k) = Vc;
    Npn_use_hist(k) = Npn_use;
    chi_v_hist(k) = chi_v;
    chi_los_hist(k) = chi_los;
    psi_hist(k) = psi;
    phi_cmd_hist(k) = phi_cmd;
    theta_cmd_hist(k) = theta_cmd;
    alpha_hist(k) = alpha;
    alpha_cmd_hist(k) = alpha_cmd;
    CL_req_hist(k) = CL_req;
    CL_hist(k) = CL;
    a_up_aero_hist(k) = a_up_aero;
    Vcmd_hist(k) = V_cmd;

    %% 命中判据
    if R_dist < 5e3
        fprintf('Target reached at t=%.2f s, miss distance=%.1f m\n', t(k), R_dist);
        inf_time = k;
        break;
    end

    %% 数值保护
    if any(~isfinite([r;v;phi;theta;psi;Ma;qbar;alpha;beta]))
        fprintf('NaN/Inf at step %d, t=%.2f s\n', k, t(k));
        inf_time = k;
        break;
    end
end

%% 截断
idx = 1:inf_time;
tt = t(idx);

%% 图1：ECEF轨迹
figure;
plot3(Rhist(idx,1)/1e3,Rhist(idx,2)/1e3,Rhist(idx,3)/1e3,'b','LineWidth',1.2); hold on; grid on; axis equal;
xlabel('X_{ECEF} (km)'); ylabel('Y_{ECEF} (km)'); zlabel('Z_{ECEF} (km)');
title('Cruise Trajectory in ECEF (3-DOF, Anti-U-Divergence)');
[xe,ye,ze]=sphere(60);
surf(Re*xe/1e3,Re*ye/1e3,Re*ze/1e3,'FaceAlpha',0.08,'EdgeColor','none','FaceColor',[0.2 0.6 1.0]);
legend('Vehicle Trajectory','Earth');

%% 图2：距离
figure;
plot(tt, Dist_hist(idx)/1000, 'm-', 'LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('Distance to Target (km)');
title('Real-time Distance between Vehicle and Target');

%% 图3：闭合速度
figure;
plot(tt, Vc_hist(idx), 'LineWidth',1.3); grid on;
xlabel('Time (s)'); ylabel('V_c (m/s)');
title('Closing Speed (V_c > 0 means closing)');
yline(0,'r--');

%% 图4：航向
figure;
plot(tt, rad2deg(chi_v_hist(idx)), 'b', 'LineWidth',1.2); hold on; grid on;
plot(tt, rad2deg(chi_los_hist(idx)), 'r--', 'LineWidth',1.2);
plot(tt, rad2deg(psi_hist(idx)), 'k-.', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Angle (deg)');
title('Heading Diagnostics');
legend('\chi_v','\chi_{LOS}','\psi','Location','best');

%% 图5：CL需求跟踪
figure;
plot(tt, CL_req_hist(idx), 'r', 'LineWidth',1.2); hold on; grid on;
plot(tt, CL_hist(idx), 'b--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('C_L');
title('Lift Tracking');
legend('C_{L,req}','C_L','Location','best');

%% 图6：向上气动力
figure;
plot(tt, a_up_aero_hist(idx), 'LineWidth',1.2); hold on; grid on;
yline(9.81,'r--','g');
xlabel('Time (s)'); ylabel('a_{up,aero} (m/s^2)');
title('Upward Aerodynamic Acceleration');

%% 图7：指令速度调度
figure;
plot(tt, Vcmd_hist(idx), 'LineWidth',1.2); hold on; grid on;
plot(tt, vecnorm(Vhist(idx,:),2,2), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Speed (m/s)');
title('Speed Command Schedule');
legend('V_{cmd}','|V|','Location','best');

%% 图8：模式切换
figure;
stairs(tt, mode_hist(idx), 'LineWidth',1.5); grid on;
xlabel('Time (s)'); ylabel('Mode');
title('Guidance Mode (0: PN-dominant, 1: PP-dominant)');
ylim([-0.2,1.2]);