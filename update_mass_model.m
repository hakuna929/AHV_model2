function [m_new, mf_new, mdot_f] = update_mass_model( ...
    m_struct, mf_old, Ma, alpha, beta, H, dT, dt, ...
    use_alpha_deg, mdot_min, mdot_max)
%UPDATE_MASS_MODEL 质量消耗模型
%
% 输入：
%   m_struct        - 结构质量 (kg)
%   mf_old          - 当前燃料质量 (kg)
%   Ma              - 马赫数
%   alpha           - 攻角（rad）
%   beta            - 侧滑角（rad）
%   H_km            - 高度（km）
%   dT              - 油门/推力控制量
%   dt              - 时间步长 (s)
%   use_alpha_deg   - true: 公式中的 alpha 以度计；false: 以弧度计
%   mdot_min        - 燃料质量流率下限 (kg/s)
%   mdot_max        - 燃料质量流率上限 (kg/s)
%
% 输出：
%   m_new           - 更新后的总质量 (kg)
%   mf_new          - 更新后的燃料质量 (kg)
%   mdot_f          - 实际燃料消耗率 (kg/s)

    if nargin < 10 || isempty(use_alpha_deg)
        use_alpha_deg = true;
    end
    if nargin < 11 || isempty(mdot_min)
        mdot_min = 0.0;
    end
    if nargin < 12 || isempty(mdot_max)
        mdot_max = inf;
    end

    if use_alpha_deg
        alpha_use = alpha * 180/pi;
    else
        alpha_use = alpha;
    end

    cb = cos(beta);

    % 质量消耗模型（按你给出的公式）
    mdot_per_dT = ...
        2.4805 ...
        - 0.05455 * alpha_use ...
        + 0.001599 * alpha_use^2 ...
        - 0.204 * H ...
        + 0.486 * Ma * cb ...
        + 0.002515 * alpha_use * H ...
        + 0.003635 * H^2 ...
        + 0.008598 * (Ma * cb) * alpha_use ...
        - 0.01216 * (Ma * cb) * H;

    mdot_f = dT * mdot_per_dT;

    if ~isfinite(mdot_f)
        mdot_f = 0.0;
    end

    mdot_f = min(max(mdot_f, mdot_min), mdot_max);

    mf_new = mf_old - mdot_f * dt;
    mf_new = max(mf_new, 0.0);

    m_new = m_struct + mf_new;

    % 由于上面可能把燃料耗尽，这里再反算一次真实消耗率
    mdot_f = (mf_old - mf_new) / max(dt,1e-6);
end