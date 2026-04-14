clear; clc;

%% 时间
dt = 0.001;
T  = 400;
t  = 0:dt:T;
N  = numel(t);

%% 地球常量
mu  = 3.986004418e14;
J2  = 1.08263e-3;
Re  = 6378137;
we  = 7.2921150e-5;
omega_ie = [0;0;we];

%% 飞行器参数
m     = 671.33;
S_ref = 0.2986;
b_ref = 0.8;
c_ref = 0.3732;

Ixz = 430;
I = [34.13 0 -Ixz;
     0    1040 0;
    -Ixz  0 1034];
Iinv = inv(I);

%% 定速定高巡航初始条件（ECI）
h0   = 30e3;        %高度
lat0 = 19.2;   %纬度海南文昌
lon0 = 110.5;  %经度海南文昌
r = lla2ecef_wgs84(lat0, lon0, h0);
r0_ecef = r;

Ma_r = 6.5; %期望巡航马赫数
[a0,rh0_0] = atmos_simple(h0);
V0 = Ma_r * a0;   %巡航速度


chi0 = 30*pi/180; %航迹方位角/航向角 
gamma0 = 0*pi/180;  %航迹倾斜角/飞行路径角
eE = [-sin(lon0); cos(lon0); 0];
eN = [-sin(lat0)*cos(lon0); -sin(lat0)*sin(lon0); cos(lat0)];
eU = [cos(lat0)*cos(lon0); cos(lat0)*sin(lon0); sin(lat0)];
v_rel0 = V0*(cos(gamma0)*cos(chi0)*eN + cos(gamma0)*sin(chi0)*eE + sin(gamma0)*eU);
v = v_rel0 + cross(omega_ie,r);

% 姿态与角速率
% phi=0; theta=0; psi=chi0; p=0; q=0; rrate=0;
% 姿态与角速率初始化 (确保 alpha=0, beta=0) 
v_air_eci_0 = v_rel0; % 因为 v = v_rel0 + cross(omega_ie,r0_ecef)
Vair_0 = norm(v_air_eci_0);
xb = v_air_eci_0 / Vair_0;

% 2. 假设初始无滚转(机翼水平)，机体系 Y 轴必须垂直于 X 轴和当地垂线(eU)
yb = cross(eU, xb);
yb = yb / norm(yb);

% 3. 根据右手法则确定机体系 Z 轴
zb = cross(xb, yb);

% 4. 组装初始的 "机体系到ECI系" 旋转矩阵 (对应你代码里的 Ci_b)
C_b2i_0 = [xb, yb, zb];

% 5. 从旋转矩阵反解出 ECI 系下的初始欧拉角 (基于3-2-1顺规)
theta = asin(-C_b2i_0(3,1));
psi   = atan2(C_b2i_0(2,1), C_b2i_0(1,1));
phi   = atan2(C_b2i_0(3,2), C_b2i_0(3,3));

% 初始化角速率
p = 0; q = 0; rrate = 0;


% 控制初值
de1=0; de2=0; dr=0; dT=1;

%% 目标点设置（与起点球面距离）
R0 = Re + h0;
s_target = 1000e3;
sigma = s_target / R0;
az0 = chi0;
n_gc = cos(az0)*eN + sin(az0)*eE;
rT = R0*(cos(sigma)*(r0_ecef/R0) + sin(sigma)*n_gc);
vT = [0;0;0];

%% 控制参数
Npn = 3.0;
a_cmd_max = 30;             % m/s^2

Kp_phi=2.0; Kd_p=1.0;
Kp_th =3.0; Kd_q=1.2;
Kp_psi=1.2; Kd_r=0.8;

de_lim = 20*pi/180;
dr_lim = 25*pi/180;

V_cmd = 2500;
Kv = 0.0015;

%% 记录
Rhist = zeros(N,3);
Vhist = zeros(N,3);
Hhist = zeros(N,1);

for k=1:N
    %% 地理量
    rn = norm(r); ur = r/rn;
    [lat, lon, h] = ecef2lla_wgs84(r');  % 注意 r 需为行向量
    Hhist(k)=h;

    %% 大气相对速度
    v_air_eci = v - cross(omega_ie,r);
    Vair = max(norm(v_air_eci),1e-3);

    %% 姿态矩阵
    cph=cos(phi); sph=sin(phi); cth=cos(theta); sth=sin(theta); cps=cos(psi); sps=sin(psi);
    Ci_b = [cth*cps, sph*sth*cps-cph*sps, cph*sth*cps+sph*sps;
            cth*sps, sph*sth*sps+cph*cps, cph*sth*sps-sph*cps;
            -sth,    sph*cth,             cph*cth];
    Cb_i = Ci_b';

    %% 空速到机体系
    v_air_b = Cb_i * v_air_eci;
    u=v_air_b(1); vv=v_air_b(2); w=v_air_b(3);
    alpha = atan2(w,u);
    beta  = asin(max(min(vv/Vair,1),-1));

    % 保护限幅（防数值发散）
    alpha = max(min(alpha, 10*pi/180), -2*pi/180);
    beta  = max(min(beta,  6*pi/180), -6*pi/180);

    %% 大气、马赫、动压
    [a,rho] = atmos_simple(max(h,0));
    Ma = Vair/max(a,1e-3);
    qbar = 0.5*rho*Vair^2;

    %% 重力（先算，供制导用）
    x=r(1); y=r(2); z=r(3); rr=max(norm(r),Re+1);
    g0 = -mu/rr^3*[x;y;z];
    gJ2 = [3*mu*J2*Re^2/(2*rr^5)*x*(1-5*(z/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*y*(1-5*(z/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*z*(3-5*(z/rr)^2)];
    g_eci = g0 + gJ2;
    gmag = max(norm(g_eci),1e-3);

    %% ========== 比例导引（PN） ==========
    r_rel = rT - r;
    R = max(norm(r_rel),1);
    u_los = r_rel / R;
    v_rel = vT - v;
    Vc = -dot(r_rel,v_rel)/R;                         % 闭合速度
    omega_los = cross(r_rel,v_rel)/(R^2 + 1e-6);
    a_cmd_eci = Npn * Vc * cross(omega_los, u_los);

    acn = norm(a_cmd_eci);
    if acn > a_cmd_max
        a_cmd_eci = a_cmd_eci/acn * a_cmd_max;
    end

    u_up = r/norm(r);
    a_vert = dot(a_cmd_eci,u_up);
    a_h_eci = a_cmd_eci - a_vert*u_up;
    a_h = norm(a_h_eci);

    eE_now = [-sin(lon); cos(lon); 0];
    eN_now = [-sin(lat)*cos(lon); -sin(lat)*sin(lon); cos(lat)];

    if a_h > 1e-6
        dir_h = a_h_eci/a_h;
        chi_cmd = atan2(dot(dir_h,eE_now), dot(dir_h,eN_now));
    else
        chi_cmd = psi;
    end

    theta_cmd = atan2(a_vert,gmag);
    theta_cmd = max(min(theta_cmd,10*pi/180),-10*pi/180);
    phi_cmd   = atan2(a_h,gmag);
    phi_cmd   = max(min(phi_cmd,35*pi/180),-35*pi/180);
    psi_cmd   = chi_cmd;

    %% ========== 姿态控制（PD） ==========
    e_psi = atan2(sin(psi_cmd-psi), cos(psi_cmd-psi));
    e_phi = phi_cmd - phi;
    e_th  = theta_cmd - theta;

    L_cmd = Kp_phi*e_phi - Kd_p*p;
    M_cmd = Kp_th *e_th  - Kd_q*q;
    N_cmd = Kp_psi*e_psi - Kd_r*rrate;

    de_sym  = M_cmd;
    de_diff = L_cmd;
    de1 = 0.5*(de_sym + de_diff);
    de2 = 0.5*(de_sym - de_diff);
    dr  = N_cmd;

    de1 = max(min(de1,de_lim),-de_lim);
    de2 = max(min(de2,de_lim),-de_lim);
    dr  = max(min(dr,dr_lim),-dr_lim);

    % 推力速度环
    dT = dT + Kv*(V_cmd - Vair);
    dT = max(min(dT,1.0),0.05);

    %% 气动与推力
    [CL,CD,CY,Cl,Cm,Cn] = aero_coeffs(Ma,alpha,beta,de1,de2,dr,dT);
    L = CL*qbar*S_ref; D=CD*qbar*S_ref; Y=CY*qbar*S_ref;
    Fw = [-D;Y;-L];

    ca=cos(alpha); sa=sin(alpha); cb=cos(beta); sb=sin(beta);
    Cb_w = [ ca*cb, -ca*sb, -sa;
             sb,     cb,     0;
             sa*cb, -sa*sb,  ca ];
    Fa_b = Cb_w*Fw;

    CT = thrust_coeffs(Ma,alpha,beta,dT);
    T_eng = CT*qbar*S_ref;
    Ft_b = [T_eng;0;0];

    Mb = [Cl*qbar*S_ref*b_ref;
          Cm*qbar*S_ref*c_ref;
          Cn*qbar*S_ref*b_ref];

    %% 平动
    F_eci = Ci_b*(Fa_b + Ft_b);
    a_eci = g_eci + F_eci/m;
    v = v + a_eci*dt;
    r = r + v*dt;

    %% 转动
    omega_b = [p;q;rrate];
    omegadot = Iinv*(Mb - cross(omega_b,I*omega_b));
    p = p + omegadot(1)*dt;
    q = q + omegadot(2)*dt;
    rrate = rrate + omegadot(3)*dt;

    cth_safe = cos(theta);
    if abs(cth_safe) < 1e-3
        cth_safe = sign(cth_safe)*1e-3;
        if cth_safe == 0, cth_safe = 1e-3; end
    end
    E = [1, sph*tan(theta), cph*tan(theta);
         0, cph,           -sph;
         0, sph/cth_safe,   cph/cth_safe];
    eulerdot = E*omega_b;
    phi = phi + eulerdot(1)*dt;
    theta = theta + eulerdot(2)*dt;
    psi = psi + eulerdot(3)*dt;

    %% 记录
    Rhist(k,:)=r.';
    Vhist(k,:)=v.';

    %% 命中判据
    if R < 5e3
        fprintf('Target reached at t=%.2f s, miss distance=%.1f m\n', t(k), R);
    
        Rhist = Rhist(1:k,:);
        Vhist = Vhist(1:k,:);
        Hhist = Hhist(1:k,:);
        break;
    end

    %% NaN保护
    if any(~isfinite([r;v;phi;theta;psi;p;q;rrate;Ma;qbar;alpha;beta]))
        fprintf('NaN/Inf at step %d, t=%.2f s\n',k,t(k));
        inf_time = k ;
        break;
    end
end

%% ECI轨迹
figure;
plot3(Rhist(1:inf_time,1)/1e3,Rhist(1:inf_time,2)/1e3,Rhist(1:inf_time,3)/1e3,'b','LineWidth',1.2); hold on; grid on; axis equal;
xlabel('X_{ECI} (km)'); ylabel('Y_{ECI} (km)'); zlabel('Z_{ECI} (km)');
title('Cruise Trajectory in ECI');

[xe,ye,ze]=sphere(60);
surf(Re*xe/1e3,Re*ye/1e3,Re*ze/1e3,'FaceAlpha',0.08,'EdgeColor','none','FaceColor',[0.2 0.6 1.0]);
legend('Vehicle Trajectory','Earth');