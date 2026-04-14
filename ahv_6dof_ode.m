function dX = ahv_6dof_ode(t, X, U, p)
% X: [x y z u v w phi theta psi p q r mf]^T
% U: [de1 de2 dr dT df]^T

% 解包状态
x    = X(1);  y = X(2);  z = X(3);
u    = X(4);  v = X(5);  w = X(6);
phi  = X(7);  theta = X(8);  psi = X(9);
p_r  = X(10); q = X(11); r = X(12);
mf   = X(13);

% 解包输入
de1 = U(1);
de2 = U(2);
dr  = U(3);
dT  = U(4);
df  = U(5); 

% 当前质量
m = p.m + mf;

% 速度与姿态
V   = sqrt(u^2+v^2+w^2);
alpha = atan2(w,u);
beta  = atan2(v, sqrt(u^2+w^2));

% 高度（NED）
H = -z;

% 大气
[a, rho] = atmos_simple(H);
Ma = V / max(a,1);   % 避免除零

% 动压
qbar = 0.5*rho*V^2;

% 计算气动力系数（根据你给的多项式）
[CL, CD, CY, Cl, Cm, Cn] = aero_coeffs(Ma, alpha, beta, de1, de2, dr, dT, H);

% 计算推力和燃油流量
[Tx_body, Tz_body, mdot_f] = propulsion_model(Ma, alpha, beta, H, dT, p);

% 气动力（机体坐标）
L = CL * qbar * p.S;
D = CD * qbar * p.S;
Yf = CY * qbar * p.S;

% 机体轴向/侧向/法向力
% 约定：X 轴向前, Z 向下，升力向 -Z
Fx_aero = -D*cos(alpha) + -L*sin(alpha);   % 近似
Fz_aero = -D*sin(alpha) +  L*cos(alpha);
Fy_aero =  Yf;

% 推力（已经在机体坐标）
Fx = Fx_aero + Tx_body;
Fy = Fy_aero;
Fz = Fz_aero + Tz_body;

% 重力在机体坐标
g = 9.80665;
Rnb = angle2dcm(psi, theta, phi, 'ZYX');   % from body to NED
Rbn = Rnb.';                               % NED->body
Fg_ned = [0; 0; m*g];                      % NED: z 向下
Fg_b   = Rbn * Fg_ned;

Fx = Fx + Fg_b(1);
Fy = Fy + Fg_b(2);
Fz = Fz + Fg_b(3);

% 力矩（气动力矩）
[Cl_b, Cm_b, Cn_b] = deal(Cl, Cm, Cn);   % 直接复用 aero_coeffs 的输出，避免缺失函数
% 这里直接用给定的 C_l, C_m, C_n 系数
L_m = Cl * qbar * p.S * p.b;
M_m = Cm * qbar * p.S * p.c;
N_m = Cn * qbar * p.S * p.b;

% 刚体转动方程
Ix  = p.Ix; Iy = p.Iy; Iz = p.Iz; Ixz = p.Ixz;

% 考虑 Ixz 的耦合项
Gamma = Ix*Iz - Ixz^2;
Gamma1 = (Ixz*(Ix-Iy+Iz))/Gamma;
Gamma2 = (Iz*(Iz-Iy)+Ixz^2)/Gamma;
Gamma3 = Ix/Gamma;
Gamma4 = Ixz/Gamma;
Gamma5 = (Iz-Ix)/Iy;
Gamma6 = Ixz/Iy;
Gamma7 = ((Ix-Iy)*Ix+Ixz^2)/Gamma;
Gamma8 = Ix/Gamma;

lp = L_m; mp = M_m; np = N_m;

p_dot = Gamma1*p_r*q + Gamma2*q*r + Gamma3*lp + Gamma4*np;
q_dot = Gamma5*p_r*r - Gamma6*(p_r^2 - r^2) + mp/Iy;
r_dot = Gamma7*p_r*q + Gamma1*q*r + Gamma4*lp + Gamma8*np;

% 线加速度
u_dot = r*v - q*w + Fx/m;
v_dot = p_r*w - r*u + Fy/m;
w_dot = q*u - p_r*v + Fz/m;

% 位置微分：速度从机体坐标转换到 NED
Vb = [u; v; w];
Vn = Rnb * Vb;
x_dot = Vn(1);
y_dot = Vn(2);
z_dot = Vn(3);

% 姿态角微分
T = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
     0 cos(phi)           -sin(phi);
     0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
euler_dot = T * [p_r; q; r];
phi_dot   = euler_dot(1);
theta_dot = euler_dot(2);
psi_dot   = euler_dot(3);

% 燃油质量变化
mf_dot = -mdot_f;

dX = [x_dot; y_dot; z_dot;
      u_dot; v_dot; w_dot;
      phi_dot; theta_dot; psi_dot;
      p_dot; q_dot; r_dot;
      mf_dot];
end