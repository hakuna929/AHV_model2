%% flight_corridor_dQ_qbar.m
% 飞行走廊：动压约束 + 热流约束
% 需求：网格阴影填充 + 3条约束线（qmin/qmax/dQmax），不画 Vmin/Vmax
clear; clc; close all;

%% ===================== 基本参数 =====================
% 速度、高度范围
V_vec = (1.0e3:10:2.5e3)';      % m/s
h_vec = (25e3:0.01e3:40e3)';        % m

% 动压约束
qbar_min = 2e4;      % Pa
qbar_max = 9e4;      % Pa

% 热流约束参数（按图中公式）
C1   = 13500;
n    = 0.5;
mexp = 3.15;
Rd   = 0.12;         % m
rho0 = 1.225;        % kg/m^3
g0   = 9.8066;       % m/s^2
R0   = 6371393;      % m

% 热流上限
dQ_max = 7.9e4;      % W/m^2

%% ===================== 存储矩阵 =====================
qbar_map = nan(numel(h_vec), numel(V_vec));
dQ_map   = nan(numel(h_vec), numel(V_vec));
feasible = false(numel(h_vec), numel(V_vec));

%% ===================== 计算走廊 =====================
for ih = 1:numel(h_vec)
    h = h_vec(ih);

    % 大气参数
    [~, rho] = atmos_simple(h);

    for iv = 1:numel(V_vec)
        V = V_vec(iv);

        % 动压
        qbar = 0.5 * rho * V^2;
        qbar_map(ih, iv) = qbar;

        % 热流 dQ：按给定公式
        dQ = (C1 / sqrt(Rd)) ...
            * (rho / rho0)^n ...
            * (V / sqrt(g0 * R0))^mexp;
        dQ_map(ih, iv) = dQ;

        % 可行域判断
        feasible(ih, iv) = (qbar >= qbar_min) && (qbar <= qbar_max) && (dQ <= dQ_max);
    end
end

%% ===================== 绘图：飞行走廊 =====================
figure('Color','w','Name','Flight Corridor'); hold on; box on; grid on;

% ---------- 网格阴影填充（只对 feasible 区域着色） ----------
% 这里用两色 colormap：0->白色，1->浅黄色
hImg = imagesc(V_vec/1e3, h_vec/1e3, double(feasible));
set(gca, 'YDir', 'normal');

% 颜色：不可行(0)=白色，可行(1)=浅黄色
colormap([1 1 1; 0.95 0.9 0.65]);

% 关键：只让 feasible=1 的网格有透明度；feasible=0 完全透明
shadeAlpha = 0.55;                 % 阴影透明度（0~1）
set(hImg, 'AlphaData', shadeAlpha * double(feasible));

% ---------- 3条约束线 ----------
% q = qbar_min
contour(V_vec/1e3, h_vec/1e3, qbar_map, [qbar_min qbar_min], ...
    'b', 'LineWidth', 2);

% q = qbar_max
contour(V_vec/1e3, h_vec/1e3, qbar_map, [qbar_max qbar_max], ...
    'r', 'LineWidth', 2);

% dQ = dQ_max
contour(V_vec/1e3, h_vec/1e3, dQ_map, [dQ_max dQ_max], ...
    'k--', 'LineWidth', 2);

% 坐标和标题
xlabel('Velocity / (km s^{-1})', 'FontSize', 12);
ylabel('Altitude / km', 'FontSize', 12);
title('Flight Corridor of AHV', 'FontSize', 13);
xlim([1.5 3.0]);
ylim([0 40]);

% ---------- 文本标注（可选） ----------
text(1.75, 35, sprintf('q = %.0f \\times 10^4 Pa', qbar_min/1e4), ...
    'FontSize', 11, 'Color', 'b');
text(1.75, 17, sprintf('q = %.0f \\times 10^4 Pa', qbar_max/1e4), ...
    'FontSize', 11, 'Color', 'r');
text(2.05, 12, sprintf('dQ = %.1e W/m^2', dQ_max), ...
    'FontSize', 11, 'Color', 'k');

% ---------- 图例（用虚拟对象） ----------
hLeg1 = plot(nan, nan, 'b-',  'LineWidth', 2);
hLeg2 = plot(nan, nan, 'r-',  'LineWidth', 2);
hLeg3 = plot(nan, nan, 'k--', 'LineWidth', 2);
hLeg4 = patch(nan, nan, [0.95 0.9 0.65], 'FaceAlpha', shadeAlpha, 'EdgeColor', 'none');

legend([hLeg1, hLeg2, hLeg3, hLeg4], ...
    {'q_{min} dynamic pressure limit', ...
     'q_{max} dynamic pressure limit', ...
     'Heating rate limit', ...
     'Feasible flight corridor'}, ...
    'Location', 'northwest');

view(2);