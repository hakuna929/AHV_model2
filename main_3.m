% 最短时间轨迹优化（direct transcription + fmincon）
% 状态: x=[r;v] (ECEF, 6x1)
% 控制: u=[alpha;phi;dT] (rad, rad, throttle)
% 约束: 高度、动压、迎角；终端距离

clear; clc;

%% ====== 常量/参数（沿用你 main_2.m 的设定）======
mu  = 3.986004418e14;
J2  = 1.08263e-3;
Re  = 6378137.0;
we  = 7.2921150e-5;
omega_ie = [0;0;we];

m     = 671.33;
S_ref = 0.2986;

% 舵面固定
de1=0; de2=0; dr=0;

% 初始点
h0   = 30e3;
lat0 = 19.2;
lon0 = 110.5;
r0 = lla2ecef_cgcs2000(lat0, lon0, h0); r0 = r0(:);

Ma0 = 6.5;
[a0,~] = atmos_simple(h0);
V0 = Ma0*a0;

% 目标点
lat_T = 13.35;
lon_T = 144.55;
h_T   = 30e3;
rT = lla2ecef_cgcs2000(lat_T, lon_T, h_T); rT = rT(:);

% 初始速度方向（沿用你原逻辑：大圆切线方向）
r0_hat = r0 / norm(r0);
v_dir = rT - (rT' * r0_hat) * r0_hat;
v_dir = v_dir / norm(v_dir);
v0 = V0 * v_dir;

x0 = [r0; v0];

%% ====== 优化设置 ======
K = 301;                 % 节点数（先别太大）
Rtol = 5e3;              % 终端距离阈值 (m)

% 约束（你可按需求改）
h_min = 25e3;  h_max = 35e3;      % 高度走廊
q_max = 60e3;                     % 动压上限 (Pa) 需要你按模型合理设
alpha_min = -2*pi/180;
alpha_max = 15*pi/180;
phi_min = -60*pi/180;
phi_max =  60*pi/180;
dT_min = 0.1;  dT_max = 1.0;

tf_min = 10;          % 防止 tf->0
tf_max = 4000;

%% ====== 初值（非常关键：建议用你现有PN仿真结果做初值）
% 这里给一个"直线插值"的弱初值，能跑但收敛可能差
tf0 = 1500;

X0 = zeros(6,K);
for k=1:K
    tau = (k-1)/(K-1);
    X0(:,k) = [ (1-tau)*r0 + tau*rT; v0 ]; % 速度先保持v0
end

U0 = zeros(3,K);
U0(1,:) = 2*pi/180;     % alpha 初值
U0(2,:) = 0;            % phi 初值
U0(3,:) = 1.0;          % dT 初值

z0 = packZ(X0,U0,tf0);

%% ====== 变量边界 ======
lb = -inf(size(z0)); ub = inf(size(z0));

% X不做硬边界（可选加速度/速度界限也能加）
% U边界
[idx_alpha, idx_phi, idx_dT, idx_tf] = indexMap(K);

lb(idx_alpha) = alpha_min;  ub(idx_alpha) = alpha_max;
lb(idx_phi)   = phi_min;    ub(idx_phi)   = phi_max;
lb(idx_dT)    = dT_min;     ub(idx_dT)    = dT_max;

lb(idx_tf) = tf_min; ub(idx_tf) = tf_max;

%% ====== fmincon ======
opt = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations',2e6, ...
    'MaxIterations',200, ...
    'Display','iter', ...
    'ConstraintTolerance',1e-3, ...
    'OptimalityTolerance',1e-3, ...
    'StepTolerance',1e-6);

problem.objective = @(z) objective_tf(z,K);
problem.x0 = z0;
problem.lb = lb;
problem.ub = ub;
problem.nonlcon = @(z) nonlcon_traj(z,K,x0,rT,Rtol, ...
    h_min,h_max,q_max,alpha_min,alpha_max, ...
    mu,J2,Re,omega_ie,m,S_ref,de1,de2,dr);
problem.solver = 'fmincon';
problem.options = opt;

[zsol, fval, exitflag, output] = fmincon(problem);

%% ====== 解包&绘图 ======
[Xsol,Usol,tfsol] = unpackZ(zsol,K);
fprintf('Solved tf = %.2f s\n', tfsol);

r_sol = Xsol(1:3,:);
v_sol = Xsol(4:6,:);

figure; plot3(r_sol(1,:)/1e3,r_sol(2,:)/1e3,r_sol(3,:)/1e3,'b-'); grid on; axis equal;
xlabel('X km'); ylabel('Y km'); zlabel('Z km'); title('Optimized trajectory (ECEF)');

% 距离曲线
dist = vecnorm(r_sol - rT,2,1);
tgrid = linspace(0,tfsol,K);
figure; plot(tgrid, dist/1e3,'m-','LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('dist to target (km)'); title('Distance');

% alpha/phi/dT
figure;
subplot(3,1,1); plot(tgrid, Usol(1,:)*180/pi); grid on; ylabel('\alpha (deg)');
subplot(3,1,2); plot(tgrid, Usol(2,:)*180/pi); grid on; ylabel('\phi (deg)');
subplot(3,1,3); plot(tgrid, Usol(3,:)); grid on; ylabel('dT'); xlabel('t (s)');

%% ====== helper: pack/unpack/index/obj/constraints ======
function z = packZ(X,U,tf)
    z = [X(:); U(:); tf];
end

function [X,U,tf] = unpackZ(z,K)
    nx=6; nu=3;
    X = reshape(z(1:nx*K), nx, K);
    U = reshape(z(nx*K+1 : nx*K+nu*K), nu, K);
    tf = z(end);
end

function J = objective_tf(z,K)
    tf = z(end);
    J = tf;
end

function [idx_alpha, idx_phi, idx_dT, idx_tf] = indexMap(K)
    nx=6; nu=3;
    baseX = 0;
    baseU = nx*K;
    idx_alpha = baseU + (1:nu:nu*K);
    idx_phi   = baseU + (2:nu:nu*K);
    idx_dT    = baseU + (3:nu:nu*K);
    idx_tf    = nx*K + nu*K + 1;
end

function [c,ceq] = nonlcon_traj(z,K,x0,rT,Rtol, ...
    h_min,h_max,q_max,alpha_min,alpha_max, ...
    mu,J2,Re,omega_ie,m,S_ref,de1,de2,dr)

    [X,U,tf] = unpackZ(z,K);
    dt = tf/(K-1);

    % 初值固定：X(:,1)=x0
    ceq = X(:,1) - x0;

    % 动力学离散约束
    for k=1:K-1
        xk = X(:,k);
        uk = U(:,k);
        xkp1 = X(:,k+1);

        xdot = dyn_ecef_3dof(xk,uk,mu,J2,Re,omega_ie,m,S_ref,de1,de2,dr);
        ceq = [ceq; xkp1 - xk - dt*xdot];
    end

    % 终端距离约束: ||rK-rT|| <= Rtol
    rK = X(1:3,end);
    c_term = norm(rK - rT) - Rtol;

    % 路径约束（高度、动压）
    c_path = [];
    for k=1:K
        rk = X(1:3,k);
        vk = X(4:6,k);
        alpha = U(1,k); %#ok<NASGU>

        % 高度
        [~,~,h] = ecef2lla_cgcs2000(rk');  % 注意你函数输入是行向量
        % 动压
        Vair = max(norm(vk),1e-3);
        [a,rho] = atmos_simple(max(h,0));
        qbar = 0.5*rho*Vair^2;

        c_path = [c_path;
                  h_min - h;   % <=0 代表 h>=h_min
                  h - h_max;   % <=0 代表 h<=h_max
                  qbar - q_max]; % <=0
    end

    % 迎角边界已经放在 lb/ub，这里不再写（避免重复）
    c = [c_term; c_path];
end

function xdot = dyn_ecef_3dof(x,u,mu,J2,Re,omega_ie,m,S_ref,de1,de2,dr)
    % 状态/控制
    r = x(1:3);
    v = x(4:6);
    alpha = u(1);
    phi   = u(2);
    dT    = u(3);

    % 地理量
    [lat, lon, h] = ecef2lla_cgcs2000(r');  % lat/lon 单位需与你函数一致（你主程序里用 sin(lat) 看起来是弧度）

    % 重力(J2)
    xE=r(1); yE=r(2); zE=r(3); rr=max(norm(r),Re+1);
    g0 = -mu/rr^3*[xE;yE;zE];
    gJ2 = [3*mu*J2*Re^2/(2*rr^5)*xE*(1-5*(zE/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*yE*(1-5*(zE/rr)^2);
           3*mu*J2*Re^2/(2*rr^5)*zE*(3-5*(zE/rr)^2)];
    g_ecef = g0 + gJ2;

    % 大气/动压/Ma
    v_air = v;
    Vair = max(norm(v_air),1e-3);
    [a,rho] = atmos_simple(max(h,0));
    Ma = Vair/max(a,1e-3);
    Ma = max(min(Ma, 7.0), 0.0);
    qbar = 0.5*rho*Vair^2;

    % 构造姿态：这里用"航向=速度方向 + 滚转=phi + 迎角=alpha"的简化
    % 你原代码是从(phi,theta,psi)构造Ce_b。这里为了最小框架，只需要能把风轴力转到ECEF。
    % 建议你后续把 main_2.m 的 Cned_b/Ce_ned/Ce_b 整段复制进来，显式用 psi/theta 由速度方向推出来。

    % --- 简化做法：把气动力直接按风轴->ECEF（需要你补全：Cb_w, Ce_b 等）---
    beta = 0; % 暂时先设0，后续可像 main_2.m 那样由速度在机体系计算

    [CL,CD,CY,~,~,~] = aero_coeffs(Ma,alpha,beta,de1,de2,dr,dT);
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

    % 这里缺 Ce_b（body->ecef）。你可以直接移植 main_2.m 里那套：
    % 1) 从(phi,theta,psi)得Cned_b
    % 2) 从(lat,lon)得Ce_ned
    % 3) Ce_b = Ce_ned*Cned_b
    % 在"最小版本"里先假设 body轴与ECEF一致（不真实，但框架能跑通）
    Ce_b = eye(3);

    F_ecef = Ce_b*(Fa_b + Ft_b);

    a_ecef = g_ecef + F_ecef/m - 2*cross(omega_ie, v) - cross(omega_ie, cross(omega_ie, r));

    xdot = [v; a_ecef];
end