function [a, rho] = atmos_simple(h)
% 标准大气模型（基于几何高度分段，严格按图片公式）
% 输入: h - 几何高度，单位 m
% 输出: a - 声速，m/s； rho - 密度，kg/m^3

%% 常数定义
R0_km = 6356.766;      % 地球半径，km
p0 = 101325;           % 海平面压强，Pa
rho0 = 1.225;          % 海平面密度，kg/m^3
R = 287;               % 空气气体常数，J/(kg·K)
gamma = 1.4;           % 比热比

%% 转换为千米
h_km = h / 1000;


% -------------------------------------------------------------------------
% 20–47 km 平滑插值（PCHIP）
% 说明：
% - 端点取原分段模型在 20 km 与 47 km 的曲线值
% - 在 (20 km, 47 km) 区间内对 log(p)、T、log(rho) 做 PCHIP 插值
% - 这样能保证 p、rho 始终为正，且整体单调/形状保持
% -------------------------------------------------------------------------
if h_km > 20 && h_km < 47
    % 计算端点（使用原始分段公式计算 p/T/rho）
    [p20, T20, rho20] = atmos_piecewise(20000, R0_km, p0, rho0);
    [p47, T47, rho47] = atmos_piecewise(47000, R0_km, p0, rho0);

    % PCHIP 插值节点（几何高度 m）
    h_nodes = [20000, 47000];

    % 对 log(p)、T、log(rho) 插值
    logp = pchip(h_nodes, log([p20, p47]), h);
    T    = pchip(h_nodes, [T20, T47], h);
    logr = pchip(h_nodes, log([rho20, rho47]), h);

    p   = exp(logp);
    rho = exp(logr);

else
    % 其他区间：沿用原分段模型（含 >47.35 km 外推）
    [p, T, rho] = atmos_piecewise(h, R0_km, p0, rho0);
end

%% 声速计算
a = sqrt(gamma * R * T);   % m/s

end


function [p, T, rho] = atmos_piecewise(h, R0_km, p0, rho0)

% 常数（分段模型内用到）
h_km = h / 1000;
H_km = h_km / (1 + h_km / R0_km);

if h_km <= 11.019
    % 第一段：0 ~ 11.019 km
    W = 1 - H_km / 44.3308;
    T = 288.15 * W;
    p = W^2.2559 * p0;
    rho = W^4.2559 * rho0;

elseif h_km <= 20.0631
    % 第二段：11.0191 ~ 20.0631 km
    W = exp((14.9647 - H_km) / 6.3416);
    T = 216.65;
    p = 0.11953 * W * p0;
    rho = 0.15898 * W * rho0;

elseif h_km <= 32.1619
    % 第三段：20.0631 ~ 32.1619 km
    W = 1 + (H_km - 24.9021) / 221.552;
    T = 211.552 * W;
    p = 0.025158 * W^(-3.4169) * p0;
    rho = 0.032722 * W^(-3.5169) * rho0;

elseif h_km <= 47.3501
    % 第四段：32.1619 ~ 47.3501 km
    W = 1 + (H_km - 39.7499) / 89.4107;
    T = 250.350 * W;
    p = 2.8338e-3 * W^(-12.2011) * p0;
    rho = 3.2618e-3 * W^(-13.2021) * rho0;

else
    % 超出模型范围（>47.35 km），简单外推：第四段末端值等温指数衰减
    H_end = 47.3501e3;   % 几何高度 m

    % 末端采用第四段公式在 47.3501 km 的值（这里用"几何高度→地势高度"的严格换算）
    h_end_km = 47.3501;
    H_end_km = h_end_km / (1 + h_end_km / R0_km);
    W_end = 1 + (H_end_km - 39.7499) / 89.4107;

    T_end = 250.350 * W_end;
    p_end = 2.8338e-3 * W_end^(-12.2011) * p0;
    rho_end = 3.2618e-3 * W_end^(-13.2021) * rho0;

    % 指数衰减
    H_scale = 7000;   % 标高，m（经验外推）
    p = p_end * exp(-(h - H_end) / H_scale);
    T = T_end;
    rho = rho_end * exp(-(h - H_end) / H_scale);
end

end