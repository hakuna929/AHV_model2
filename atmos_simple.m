function [a, rho] = atmos_simple(H)
% 简化大气模型（仅演示用，可改为标准大气）
% H: m

% 线性温度梯度假设 0~32 km
T0 = 288.15;
L  = -0.0065;
R  = 287;
g0 = 9.80665;

H = max(H,0);
if H <= 32000
    T = T0 + L*H;
    p = 101325 * (T/T0).^(-g0/(L*R));
else
    % 高于 32 km 做简单指数衰减
    T = 228.65;
    p = 8680 * exp(-(H-32000)/7000);
end

rho = p./(R*T);
a   = sqrt(1.4*R*T);
end