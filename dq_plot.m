% %% heating_constraint_curve.m
% % 单独绘制热流约束曲线 dQ = dQ_max
% clear; clc; close all;
% 
% % 高度范围（与原脚本一致）
% h_vec = (25e3:0.01e3:40e3)';   % m
% 
% % 热流约束参数（与原脚本一致）
% C1   = 13500;
% n    = 0.5;
% mexp = 3.15;
% Rd   = 0.12;         % m
% rho0 = 1.225;        % kg/m^3
% g0   = 9.8066;       % m/s^2
% R0   = 6371393;      % m
% dQ_max = 7.9e4;      % W/m^2
% 
% V_dQ = nan(size(h_vec));
% 
% for ih = 1:numel(h_vec)
%     h = h_vec(ih);
%     [~, rho] = atmos_simple(h);
% 
%     % 由 dQ = dQ_max 反解 V(h)
%     V_dQ(ih) = sqrt(g0*R0) * ( (dQ_max*sqrt(Rd)) / (C1*(rho/rho0)^n) )^(1/mexp);
% end
% 
% figure('Color','w','Name','Heating constraint only'); hold on; grid on; box on;
% plot(V_dQ/1e3, h_vec/1e3, 'k--', 'LineWidth', 2);
% xlabel('Velocity / (km s^{-1})');
% ylabel('Altitude / km');
% title('Heating constraint: dQ = dQ_{max}');
% ylim([25 40]);


%% heating_constraint_curve.m
% 单独绘制热流约束曲线（钝头体驻点：Sutton–Graves 相关式）
% 约束：q_dot <= q_dot_max
% 输出：在 h–V 平面上绘制 q_dot = q_dot_max 的边界曲线

clear; clc; close all;

%% ===================== 高度范围 =====================
% 与你原脚本一致（25–40 km）
h_vec = (0:10:40e3)';   % m

%% ===================== Sutton–Graves 参数 =====================
% 头部曲率半径/钝头半径（m）
Rn = 0.12;

% Sutton–Graves 常数 K（SI 单位常见取值之一；不同资料会略有差异）
% 形式：q_dot[W/m^2] = K * sqrt(rho/Rn) * V^3
K_SG = 1.74e-4;

% 热流上限（W/m^2）
qdot_max = 7.9e5;  % 79 kW/m^2

%% ===================== 反解边界速度 V(h) =====================
% qdot_max = K * sqrt(rho/Rn) * V^3
% => V = ( qdot_max / (K*sqrt(rho/Rn)) )^(1/3)
V_lim = nan(size(h_vec));

for ih = 1:numel(h_vec)
    h = h_vec(ih);
    rho = rho_smooth_25_40(h);

    denom = K_SG * sqrt(rho / Rn);
    if denom > 0
        V_lim(ih) = (qdot_max / denom)^(1/3);
    else
        V_lim(ih) = NaN;
    end
end

%% ===================== 绘图 =====================
figure('Color','w','Name','Heating constraint only (Sutton–Graves)');
hold on; grid on; box on;

plot(V_lim/1e3, h_vec/1e3, 'k-', 'LineWidth', 2);

xlabel('Velocity / (km s^{-1})');
ylabel('Altitude / km');
title('Heating constraint (Sutton–Graves): qdot = qdot_max', ...
      'Interpreter','none');

ylim([0 40]);
% 速度范围你可按需要调整；一般会落在几 km/s 内
xlim([1.5 3]);

legend({'\dot{q} = \dot{q}_{max}'}, 'Location', 'northwest');