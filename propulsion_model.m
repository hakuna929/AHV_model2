function [Tx, Tz, mdot_f] = propulsion_model(Ma, alpha, beta, H, dT, p)
% 推力方向假设在 X-Z 平面，主要沿 X 轴

adeg = alpha;
bdeg = beta;
H_km = H/1000;

% (3-11) 推力系数 C_Tx = C_Tx(Ma,α,β,δT) ：
CTx = -0.8297 + 0.2703*Ma*cos(bdeg) - 0.1133*adeg ...
      - 1.2315*dT + 0.01695*adeg*dT - 0.04602*(Ma*cos(beta))^2*dT;

% (3-12) m_f_dot / δT = f(Ma,α,β,H)
mf_over_dT = 2.4805 - 0.05455*adeg + 0.001599*adeg^2 - 0.204*H_km ...
             + 0.486*Ma*cos(beta) + 0.002515*adeg*H_km + 0.003635*H_km^2 ...
             - 0.008598*(Ma*cos(beta))*adeg - 0.01216*(Ma*cos(beta))*H_km;

mdot_f = max(mf_over_dT * dT, 0);  % kg/s

% 由推力系数 -> 推力（N）
qbar_ref = 0.5 * p.rho_ref * (p.V_ref^2); % 参考动压，用表 3‑2 中的中值
T_mag = CTx * qbar_ref * p.S;

Tx =  T_mag * cos(alpha);  % 简化
Tz = -T_mag * sin(alpha);
end