%% CL / CD / L-D / CT vs alpha
clear; clc; close all;

Ma_list = [5 5.5 6 6.5 7];
dT_list = [0.3 0.6 0.8];

adeg = (-2:0.05:10)';   % deg
alpha_label = '\alpha (deg)';

C = lines(numel(Ma_list));
legtxt = arrayfun(@(m) sprintf('Ma=%.1f', m), Ma_list, 'UniformOutput', false);

% ---- Store results for all dT: (alpha, Ma, dT) ----
CL_all = zeros(numel(adeg), numel(Ma_list), numel(dT_list));
CD_all = zeros(numel(adeg), numel(Ma_list), numel(dT_list));
LD_all = zeros(numel(adeg), numel(Ma_list), numel(dT_list));

for k = 1:numel(dT_list)
    dT = dT_list(k);

    for i = 1:numel(Ma_list)
        Ma = Ma_list(i);

        % CLa(Ma, alpha)
        CLa = 0.1498 ...
            - 0.02751 * Ma ...
            + 0.07235 * adeg ...
            - 0.003368 * adeg * Ma ...
            + 0.002343 * Ma.^2 ...
            + 0.001185 * adeg.^2;

        % C_Le = C_Le(Ma, alpha, delta_T)
        CLe = 0.7215 ...
            + 0.02635  * adeg ...
            + 0.1147   * Ma ...
            - 0.002795 * Ma * adeg ...
            - 0.5782   * sqrt(Ma) ...
            + 0.2894   * dT ...
            - 0.004363 * adeg * dT ...
            + 0.01083  * Ma * dT;

        % CDa(Ma, alpha)
        CDa = 0.05099 ...
            - 0.004863 * Ma ...
            + 0.002967 * adeg ...
            + 0.001364 * adeg.^2;

        % C_De = C_De(Ma, alpha, delta_T)
        CDe = 0.002339 * adeg ...
            + 0.00012182 * adeg.^2 ...
            - 0.00033126 * Ma * adeg ...
            + 0.005557 * adeg * dT;

        CL = CLa + CLe;
        CD = CDa + CDe;

        CL_all(:, i, k) = CL;
        CD_all(:, i, k) = CD;
        LD_all(:, i, k) = CL ./ CD;
    end

    %% -------- Figure: CL (per dT) --------
    figure('Name', sprintf('CL vs alpha (dT=%.1f)', dT)); hold on; grid on;
    for i = 1:numel(Ma_list)
        plot(adeg, CL_all(:,i,k), 'LineWidth', 1.8, 'Color', C(i,:));
    end
    xlabel(alpha_label); ylabel('C_L');
    title(sprintf('C_L vs \\alpha (dT=%.1f)', dT));
    legend(legtxt, 'Location','best');

    %% -------- Figure: CD (per dT) --------
    figure('Name', sprintf('CD vs alpha (dT=%.1f)', dT)); hold on; grid on;
    for i = 1:numel(Ma_list)
        plot(adeg, CD_all(:,i,k), 'LineWidth', 1.8, 'Color', C(i,:));
    end
    xlabel(alpha_label); ylabel('C_D');
    title(sprintf('C_D vs \\alpha (dT=%.1f)', dT));
    legend(legtxt, 'Location','best');

    %% -------- Figure: L/D (per dT) --------
    figure('Name', sprintf('L-D vs alpha (dT=%.1f)', dT)); hold on; grid on;
    for i = 1:numel(Ma_list)
        plot(adeg, LD_all(:,i,k), 'LineWidth', 1.8, 'Color', C(i,:));
    end
    xlabel(alpha_label); ylabel('C_L / C_D');
    title(sprintf('L/D vs \\alpha (dT=%.1f)', dT));
    legend(legtxt, 'Location','best');
end


%% -------- Figures: CT (separate for each dT) --------
for k = 1:numel(dT_list)
    dT = dT_list(k);

    CT_all = zeros(numel(adeg), numel(Ma_list));

    for i = 1:numel(Ma_list)
        Ma = Ma_list(i);

        % CTc(Ma, alpha, dT)
        CTc =  0.1029      * dT ...
            - 0.02022     * Ma * dT ...
            - 0.001757    * adeg  * dT ...
            + 0.000088233 * adeg.^2 * dT ...
            + 0.001221    * Ma.^2 .* dT;

        % CTn(Ma, alpha, dT)
        CTn =  0.03791 ...
            + 0.005176   * adeg ...
            - 0.01235    * Ma ...
            - 0.00054887 * Ma * adeg ...
            + 0.00096897 * Ma.^2 ...
            + 0.05627    * dT ...
            - 0.00077467 * (adeg * dT) ...
            + 0.002103   * Ma  * dT;

        CT_all(:,i) = CTc + CTn;
    end

    figure('Name', sprintf('CT vs alpha (dT=%.1f)', dT)); hold on; grid on;
    for i = 1:numel(Ma_list)
        plot(adeg, CT_all(:,i), 'LineWidth', 1.8, 'Color', C(i,:));
    end
    xlabel(alpha_label); ylabel('C_T');
    title(sprintf('C_T vs \\alpha (dT=%.1f)', dT));
    legend(legtxt, 'Location','best');
end


%% -------- Figure: CT vs dT at fixed alpha=0 deg --------
alpha0 = 4;                 % fixed alpha (deg)
dT_sweep = (0.6:0.01:1)';     % dT range, adjust if needed

CT_dT = zeros(numel(dT_sweep), numel(Ma_list));

for i = 1:numel(Ma_list)
    Ma = Ma_list(i);

    % CTc(Ma, alpha0, dT)
    CTc =  0.1029      .* dT_sweep ...
        - 0.02022     .* Ma .* dT_sweep ...
        - 0.001757    .* alpha0 .* dT_sweep ...
        + 0.000088233 .* alpha0.^2 .* dT_sweep ...
        + 0.001221    .* Ma.^2 .* dT_sweep;

    % CTn(Ma, alpha0, dT)
    CTn =  0.03791 ...
        + 0.005176   .* alpha0 ...
        - 0.01235    .* Ma ...
        - 0.00054887 .* Ma .* alpha0 ...
        + 0.00096897 .* Ma.^2 ...
        + 0.05627    .* dT_sweep ...
        - 0.00077467 .* (alpha0 .* dT_sweep) ...
        + 0.002103   .* Ma .* dT_sweep;

    CT_dT(:, i) = CTc + CTn;
end

figure('Name', sprintf('CT vs dT (alpha=%.1f deg)', alpha0)); hold on; grid on;
for i = 1:numel(Ma_list)
    plot(dT_sweep, CT_dT(:,i), 'LineWidth', 1.8, 'Color', C(i,:));
end
xlabel('dT'); ylabel('C_T');
title(sprintf('C_T vs dT at \\alpha = %.1f^\\circ', alpha0));
legend(legtxt, 'Location','best');