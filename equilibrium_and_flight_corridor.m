% 飞行走廊：动压约束 + 热流约束

clear; clc; close all;

%% ===================== 基本参数 =====================
% 速度、高度范围
V_vec = (1.5e3:10:3e3)';      % m/s
h_vec = (0:10:40e3)';          % m

% 动压约束
qbar_min = 2e4;      % Pa
qbar_max = 9e4;      % Pa

% 头部曲率半径/钝头半径（m）
Rn = 0.12;

% Sutton–Graves 常数 K（SI 单位常见取值）
K_SG = 1.74e-4;

% 热流上限
dQ_max = 7.9e5;      % W/m^2

%% ===================== 存储矩阵 =====================
qbar_map = nan(numel(h_vec), numel(V_vec));
dQ_map   = nan(numel(h_vec), numel(V_vec));
feasible = false(numel(h_vec), numel(V_vec));

%% ===================== 计算走廊 =====================
for ih = 1:numel(h_vec)
    h = h_vec(ih);
    [~, rho] = atmos_simple(h);   % 需要确保 atmos_simple 函数已定义

    for iv = 1:numel(V_vec)
        V = V_vec(iv);

        % 动压
        qbar = 0.5 * rho * V^2;
        qbar_map(ih, iv) = qbar;

        % 热流 dQ：Sutton–Graves 模型
        dQ = K_SG * sqrt(rho / Rn) * V^3;
        dQ_map(ih, iv) = dQ;

        % 可行域判断（标量）
        feasible(ih, iv) = (qbar >= qbar_min) && (qbar <= qbar_max) && (dQ <= dQ_max);
    end
end

%% ===================== 绘图 =====================
V_km = V_vec / 1e3;          % 速度 (km/s)
h_km = h_vec / 1e3;          % 高度 (km)

figure('Color','w','Name','Flight Corridor');
ax = axes('Box','on'); hold(ax,'on'); grid(ax,'on');

% 1) 填充可行域（淡黄色半透明）
contourf(ax, V_km, h_km, double(feasible), [0.5 1], ...
    'FaceColor', [0.95 0.9 0.65], 'FaceAlpha', 0.45, 'EdgeColor','none');

% 2) 绘制约束边界线
% 动压下限（蓝色实线）
[C1, h1] = contour(ax, V_km, h_km, qbar_map, [qbar_min qbar_min], ...
    'b-', 'LineWidth', 2.5);
% 动压上限（红色实线）
[C2, h2] = contour(ax, V_km, h_km, qbar_map, [qbar_max qbar_max], ...
    'r-', 'LineWidth', 2.5);
% 热流限制（黑色虚线）
[C3, h3] = contour(ax, V_km, h_km, dQ_map, [dQ_max dQ_max], ...
    'k--', 'LineWidth', 2.5);

% 3) 添加文字标注（放置于曲线合适位置）
% 自动获取每条等值线的坐标点，选择中部或末端添加文字
% 动压下限标注
x1 = C1(1,2:end); y1 = C1(2,2:end);
if ~isempty(x1)
    idx1 = round(length(x1)*0.7);  % 取曲线70%位置
    text(x1(idx1), y1(idx1), '  q_{min}=2×10^4 Pa', ...
        'Color','b', 'FontSize',10, 'FontWeight','bold', ...
        'VerticalAlignment','bottom');
end

% 动压上限标注
x2 = C2(1,2:end); y2 = C2(2,2:end);
if ~isempty(x2)
    idx2 = round(length(x2)*0.3);  % 取曲线30%位置
    text(x2(idx2), y2(idx2), '  q_{max}=9×10^4 Pa', ...
        'Color','r', 'FontSize',10, 'FontWeight','bold', ...
        'VerticalAlignment','bottom');
end

% 热流限制标注
x3 = C3(1,2:end); y3 = C3(2,2:end);
if ~isempty(x3)
    idx3 = round(length(x3)*0.5);  % 取曲线中间
    text(x3(idx3), y3(idx3), '  Heating rate limit dQ=7.9×10^5 W/m^2', ...
        'Color','k', 'FontSize',10, 'FontWeight','bold', ...
        'VerticalAlignment','top');
end

% 4) 坐标轴与外观
xlabel('Velocity / (km s^{-1})');
ylabel('Altitude / km');
title('Flight Corridor of AHV');
xlim([min(V_km), max(V_km)]);
ylim([min(h_km), max(h_km)]);
set(ax,'Layer','top');   % 网格线显示在最上层

% 5) 图例（可选，由于已有直接标注，可不添加）
% 若仍想保留图例，可加：
legend('飞行走廊','Location','northwest');