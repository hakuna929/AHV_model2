function [alpha_deg, info] = invert_alpha_from_CL(Ma, dT, CL_req, alpha_min, alpha_max)
% 依据 aero_coeffs(Ma, alpha_deg, dT) 稳健反解 alpha(deg)
% 输入:
%   Ma, dT, CL_req
%   alpha_min/max: 搜索范围(度)，建议与aero_coeffs限幅一致 [-2,10]
% 输出:
%   alpha_deg: 反解得到的alpha(度)，已在范围内
%   info: 诊断信息（是否饱和、是否单调、CL范围等）

    if nargin < 4, alpha_min = -2; end
    if nargin < 5, alpha_max = 10; end

    % --- 网格（分辨率可调：越密越准但慢一点） ---
    alpha_grid = linspace(alpha_min, alpha_max, 161);  % 161点通常足够
    CL_grid = zeros(size(alpha_grid));

    for i = 1:numel(alpha_grid)
        [CL_grid(i), ~, ~] = aero_coeffs(Ma, alpha_grid(i), dT);
    end

    % --- 找到可用的单调段（优先选择整体单调的情况） ---
    dCL = diff(CL_grid);
    mono_inc = all(dCL >= -1e-6);
    mono_dec = all(dCL <=  1e-6);

    info = struct();
    info.mono_inc = mono_inc;
    info.mono_dec = mono_dec;
    info.CL_min = min(CL_grid);
    info.CL_max = max(CL_grid);
    info.alpha_min = alpha_min;
    info.alpha_max = alpha_max;

    % --- CL_req 超范围：直接饱和到边界 ---
    if CL_req <= info.CL_min
        [~, idx] = min(CL_grid);
        alpha_deg = alpha_grid(idx);
        info.saturated = true;
        info.reason = "CL_req_below_min";
        return;
    elseif CL_req >= info.CL_max
        [~, idx] = max(CL_grid);
        alpha_deg = alpha_grid(idx);
        info.saturated = true;
        info.reason = "CL_req_above_max";
        return;
    end

    info.saturated = false;

    % --- 单调则直接插值（最干净） ---
    if mono_inc
        alpha_deg = interp1(CL_grid, alpha_grid, CL_req, 'linear');
        info.reason = "interp_monotonic_inc";
        return;
    elseif mono_dec
        alpha_deg = interp1(flip(CL_grid), flip(alpha_grid), CL_req, 'linear');
        info.reason = "interp_monotonic_dec";
        return;
    end

    % --- 非单调：用"最接近CL_req"的点做兜底，再做局部插值 ---
    [~, i0] = min(abs(CL_grid - CL_req));
    i1 = max(1, i0-1);
    i2 = min(numel(alpha_grid), i0+1);

    % 尽量找一个"跨过CL_req"的邻域（否则就直接取最近点）
    % 扩展搜索一个小窗口
    win = 10;
    ia = max(1, i0-win);
    ib = min(numel(alpha_grid), i0+win);

    alpha_win = alpha_grid(ia:ib);
    CL_win = CL_grid(ia:ib);

    % 在窗口内找是否存在跨越
    alpha_deg = alpha_grid(i0);
    info.reason = "nearest_point_fallback";

    for j = 1:numel(alpha_win)-1
        c1 = CL_win(j); c2 = CL_win(j+1);
        if (CL_req - c1) * (CL_req - c2) <= 0  % 跨越或命中
            % 线性插值
            a1 = alpha_win(j); a2 = alpha_win(j+1);
            if abs(c2 - c1) < 1e-12
                alpha_deg = 0.5*(a1+a2);
            else
                alpha_deg = a1 + (CL_req - c1) * (a2 - a1) / (c2 - c1);
            end
            info.reason = "local_bracket_interp";
            break;
        end
    end

    % 最终保护
    alpha_deg = min(max(alpha_deg, alpha_min), alpha_max);
end