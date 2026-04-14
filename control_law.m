function U = control_law(t, X, p, rt_n)
% 输出: U = [de1; de2; dr; dT; df]
% 说明: 最小可跑通版本：PN 生成加速度矢量 -> 俯仰/航向指向 + PD 自动驾驶
%
% X = [x y z u v w phi theta psi p q r mf]^T

% 状态解包
x = X(1); y = X(2); z = X(3);
u = X(4); v = X(5); w = X(6);
phi = X(7); theta = X(8); psi = X(9);
p_r = X(10); q = X(11); r = X(12);
mf = X(13);

m = p.m + mf;
g = 9.80665;

% 机体系速度 -> NED速度
Rnb = angle2dcm(psi, theta, phi, 'ZYX'); % body -> NED
Vb  = [u; v; w];
vm_n = Rnb * Vb;

% 目标（静止）
rm_n = [x; y; z];
vt_n = [0;0;0];

R = rt_n - rm_n;
V = vt_n - vm_n;

Rnorm = norm(R);
if Rnorm < 1
    Rnorm = 1;
end
rhat = R / Rnorm;

% 闭合速度
Vc = -(dot(R, V) / Rnorm);

% LOS角速度
omega_los = cross(R, V) / (Rnorm^2);

% PN参数
N = 4;              % 导航系数(可调)
a_max = 10*g;       % 最大法向加速度(可调)

% PN加速度指令 (NED)
a_cmd = N * max(Vc,0) * cross(omega_los, rhat);

% 限幅
a_norm = norm(a_cmd);
if a_norm > a_max
    a_cmd = a_cmd * (a_max / a_norm);
end

% --- 航向指向（简单版）---
dx = R(1); dy = R(2); dz = R(3);
psi_cmd = atan2(dy, dx);
psi_err = wrapToPi(psi_cmd - psi);

Kpsi = 0.5;     % 很保守，避免侧向发散（可调）
Kr   = 0.2;     % 偏航阻尼（可调）
dr = sat(Kpsi*psi_err - Kr*r, deg2rad(20));

% --- 俯仰指向 + PN上向加速度修正 ---
% 目标视线俯仰角（以水平面为基准，向上为正）
Rhor = sqrt(dx^2 + dy^2);
theta_los = atan2(-dz, max(Rhor,1));  % NED里z向下，所以 -dz 表示向上

a_up_cmd = -a_cmd(3); % NED 的 up = -z
Kpn = 0.05;           % 把PN加速度"揉"进俯仰指令的比例（可调）
theta_cmd = theta_los + Kpn*(a_up_cmd/g);

Ktheta = 2.0;         % 俯仰角P（可调）
Kq     = 0.8;         % 俯仰角速率阻尼（可调）
de = sat(Ktheta*(theta_cmd - theta) - Kq*q, deg2rad(25));

% 推力与其他舵面（先恒定）
dT = 0.2;
df = 0;

% 输出：左右升降舵对称
de1 = de;
de2 = de;

U = [de1; de2; dr; dT; df];
end

function y = sat(u, umax)
y = min(max(u, -umax), umax);
end