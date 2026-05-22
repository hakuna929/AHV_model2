% plot_mdot_f.m
% 绘制燃料秒耗量随油门开度和高度的变化
% 基于 update_mass_model.m 中的质量消耗模型
% 说明：
%   - 攻角 alpha 按"度"使用
%   - 侧滑角 beta 仍按弧度输入
%   - 高度 H 按 km 输入

clear; clc; close all;

%% =========================
%  公共参数设置
% ==========================
Ma = 6.0;                 % 马赫数
alpha_deg = 4;            % 攻角（度）
beta = 0;        % 侧滑角（弧度）
use_alpha_deg = true;     % 公式中 alpha 按"度"使用

% 注意：根据原始函数注释，H 的单位是 km
H0 = 25;                  % 固定高度（km）
dT0 = 0.5;                % 固定油门开度（0~1）

%% =========================
%  1) 燃料秒耗量随油门开度变化
% ==========================
dT_vec = linspace(0, 1, 200);
mdot_dT = zeros(size(dT_vec));

for i = 1:length(dT_vec)
    dT = dT_vec(i);
    mdot_dT(i) = calc_mdot_f(dT, Ma, alpha_deg, beta, H0, use_alpha_deg);
end

figure('Name', 'mdot_f vs dT');
plot(dT_vec, mdot_dT, 'LineWidth', 2);
grid on;
xlabel('油门开度 dT', 'Interpreter', 'tex');
ylabel('燃料秒耗量 \\dot{m}_f (kg/s)', 'Interpreter', 'tex');
title(sprintf('燃料秒耗量随油门开度变化 (H = %.2f km, Ma = %.2f, \\alpha = %.1f^\\circ)', ...
    H0, Ma, alpha_deg), 'Interpreter', 'tex');

%% =========================
%  2) 燃料秒耗量随高度变化
% ==========================
H_vec = linspace(25, 35, 200);   % 高度范围（km）
mdot_H = zeros(size(H_vec));

for i = 1:length(H_vec)
    H = H_vec(i);
    mdot_H(i) = calc_mdot_f(dT0, Ma, alpha_deg, beta, H, use_alpha_deg);
end

figure('Name', 'mdot_f vs H');
plot(H_vec, mdot_H, 'LineWidth', 2);
grid on;
xlabel('高度 H (km)', 'Interpreter', 'tex');
ylabel('燃料秒耗量 \\dot{m}_f (kg/s)', 'Interpreter', 'tex');
title(sprintf('燃料秒耗量随高度变化 (dT = %.2f, Ma = %.2f, \\alpha = %.1f^\\circ)', ...
    dT0, Ma, alpha_deg), 'Interpreter', 'tex');

%% =========================
%  3) 二维曲面：mdot_f(dT, H)
% ==========================
[dT_grid, H_grid] = meshgrid(linspace(0, 1, 80), linspace(0, 20, 80));
mdot_grid = zeros(size(dT_grid));

for i = 1:size(dT_grid, 1)
    for j = 1:size(dT_grid, 2)
        mdot_grid(i, j) = calc_mdot_f(dT_grid(i, j), Ma, alpha_deg, beta, H_grid(i, j), use_alpha_deg);
    end
end

figure('Name', 'mdot_f surface');
surf(dT_grid, H_grid, mdot_grid, 'EdgeColor', 'none');
xlabel('油门开度 dT', 'Interpreter', 'tex');
ylabel('高度 H (km)', 'Interpreter', 'tex');
zlabel('燃料秒耗量 \\dot{m}_f (kg/s)', 'Interpreter', 'tex');
title(sprintf('燃料秒耗量曲面 (Ma = %.2f, \\alpha = %.1f^\\circ)', Ma, alpha_deg), ...
    'Interpreter', 'tex');
colorbar;
grid on;
view(45, 30);

%% =========================
%  本脚本用到的局部函数
% ==========================
function mdot_f = calc_mdot_f(dT, Ma, alpha, beta, H, use_alpha_deg)
    % 计算燃料秒耗量 mdot_f
    %
    % 输入：
    %   dT              - 油门开度
    %   Ma              - 马赫数
    %   alpha           - 攻角（弧度）
    %   beta            - 侧滑角（弧度）
    %   H               - 高度（km）
    %   use_alpha_deg   - true: 公式中 alpha 以度计；false: 以弧度计
    %
    % 输出：
    %   mdot_f          - 燃料秒耗量 (kg/s)

    if use_alpha_deg
        alpha_use = alpha;
    else
         alpha_use = alpha * 180/pi;
    end

    cb = cos(beta);

    mdot_per_dT = ...
        2.4805 ...
        - 0.05455 * alpha_use ...
        + 0.001599 * alpha_use^2 ...
        - 0.204 * H ...
        + 0.486 * Ma * cb ...
        + 0.002515 * alpha_use * H ...
        + 0.003635 * H^2 ...
        - 0.008598 * (Ma * cb) * alpha_use ...
        - 0.01216 * (Ma * cb) * H;

    mdot_f = dT * mdot_per_dT;
end