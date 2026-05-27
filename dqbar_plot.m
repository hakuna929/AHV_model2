%% qbar_constraint_curves.m
% 单独绘制动压约束曲线 q = qmin 与 q = qmax
clear; clc; close all;

% 高度范围（与原脚本一致）
h_vec = (0:0.01e3:40e3)';   % m

% 动压约束（与原脚本一致）
qbar_min = 2e4;  % Pa
qbar_max = 9e4;  % Pa

V_qmin = nan(size(h_vec));
V_qmax = nan(size(h_vec));

for ih = 1:numel(h_vec)
    h = h_vec(ih);
    [~, rho] = atmos_simple(h);

    V_qmin(ih) = sqrt(2*qbar_min / rho);
    V_qmax(ih) = sqrt(2*qbar_max / rho);
end

figure('Color','w','Name','Dynamic pressure constraints only'); hold on; grid on; box on;
plot(V_qmin/1e3, h_vec/1e3, 'b-', 'LineWidth', 2);
plot(V_qmax/1e3, h_vec/1e3, 'r-', 'LineWidth', 2);
xlabel('Velocity / (km s^{-1})');
ylabel('Altitude / km');
title('Dynamic pressure constraints: q = q_{min}, q = q_{max}');
legend('q_{min}', 'q_{max}', 'Location', 'northwest');
ylim([0 40]);
% 速度范围你可按需要调整；一般会落在几 km/s 内
xlim([1.5 3]);