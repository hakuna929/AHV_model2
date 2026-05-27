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


%% -------- Figure: (CT - CD) vs alpha --------
for k = 1:numel(dT_list)
    dT = dT_list(k);

    figure('Name', sprintf('(CT-CD) vs alpha (dT=%.1f)', dT)); hold on; grid on;
    for i = 1:numel(Ma_list)
        Ma = Ma_list(i);

        % --- CD ---
        CDa = 0.05099 ...
            - 0.004863 * Ma ...
            + 0.002967 * adeg ...
            + 0.001364 * adeg.^2;

        CDe = 0.002339 * adeg ...
            + 0.00012182 * adeg.^2 ...
            - 0.00033126 * Ma * adeg ...
            + 0.005557 * adeg * dT;

        CD = CDa + CDe;

        % --- CT ---
        CTc =  0.1029      * dT ...
            - 0.02022     * Ma * dT ...
            - 0.001757    * adeg  * dT ...
            + 0.000088233 * adeg.^2 * dT ...
            + 0.001221    * Ma.^2 * dT;

        CTn =  0.03791 ...
            + 0.005176   * adeg ...
            - 0.01235    * Ma ...
            - 0.00054887 * Ma * adeg ...
            + 0.00096897 * Ma.^2 ...
            + 0.05627    * dT ...
            - 0.00077467 * (adeg * dT) ...
            + 0.002103   * Ma  * dT;

        CT = CTc + CTn;

        plot(adeg, CT - CD, 'LineWidth', 1.8, 'Color', C(i,:));
    end
    xlabel(alpha_label);
    ylabel('C_T - C_D');
    title(sprintf('C_T - C_D vs \\alpha (dT=%.1f)', dT));
    legend(legtxt, 'Location','best');
end


%% ========================================================================
%  NEW: Global max (L/D) over all alpha and dT, plus L/W and T/W vs altitude
%  NOTE:
%  - L/W = q * S * CL / (m*g)
%  - T/W = q * S * CT / (m*g)
%  You MUST set mass m and reference area S for your vehicle.
%  The ISA model below is a simple 1976-like approximation up to ~86 km.
% ========================================================================

% --- Vehicle parameters (EDIT THESE) ---
S_ref = 1.0;         % reference area [m^2]  <<< CHANGE to your vehicle
m     = 1000;        % mass [kg]             <<< CHANGE to your vehicle
g0    = 9.80665;     % gravity [m/s^2]

W = m * g0;

% --- Altitude sweep ---
h = (0:1000:50000)';     % altitude [m] (0~50 km). Adjust as needed.

% --- Find global maximum L/D for each Mach (over all alpha and dT) ---
LD_max_global    = zeros(numel(Ma_list),1);
alpha_star_g     = zeros(numel(Ma_list),1);
dT_star_g        = zeros(numel(Ma_list),1);
CL_at_star_g     = zeros(numel(Ma_list),1);
CT_at_star_g     = zeros(numel(Ma_list),1);

for i = 1:numel(Ma_list)
    LD_slice = squeeze(LD_all(:, i, :));            % [nAlpha, nDT]
    [LD_max_global(i), linIdx] = max(LD_slice(:)); % global max

    [idxAlpha, idxdT] = ind2sub(size(LD_slice), linIdx);

    alpha_star_g(i) = adeg(idxAlpha);
    dT_star_g(i)    = dT_list(idxdT);

    CL_at_star_g(i) = CL_all(idxAlpha, i, idxdT);

    % CT is not stored across dT in arrays above, so recompute at (Ma, alpha*, dT*)
    Ma = Ma_list(i);
    a  = alpha_star_g(i);
    dT = dT_star_g(i);

    CTc =  0.1029      * dT ...
        - 0.02022     * Ma * dT ...
        - 0.001757    * a  * dT ...
        + 0.000088233 * a.^2 * dT ...
        + 0.001221    * Ma.^2 .* dT;

    CTn =  0.03791 ...
        + 0.005176   * a ...
        - 0.01235    * Ma ...
        - 0.00054887 * Ma * a ...
        + 0.00096897 * Ma.^2 ...
        + 0.05627    * dT ...
        - 0.00077467 * (a * dT) ...
        + 0.002103   * Ma  * dT;

    CT_at_star_g(i) = CTc + CTn;
end

% --- Atmosphere along altitude ---
[T, a_sound, p, rho] = isa_atm_simple(h);

% Speed and dynamic pressure for each Mach at each altitude
% V(iMach, jAlt) => but we'll compute q as [nAlt, nMach]
q = zeros(numel(h), numel(Ma_list));
for i = 1:numel(Ma_list)
    V = Ma_list(i) .* a_sound;     % [m/s], vector vs h
    q(:,i) = 0.5 .* rho .* V.^2;   % dynamic pressure [Pa]
end

% --- 1) L/W vs altitude at alpha* (max L/D attack angle) ---
LW_vs_h = zeros(numel(h), numel(Ma_list));
for i = 1:numel(Ma_list)
    LW_vs_h(:,i) = (q(:,i) .* S_ref .* CL_at_star_g(i)) ./ W;
end

figure('Name','L/W vs altitude at alpha* (global max L/D)'); hold on; grid on;
for i = 1:numel(Ma_list)
    plot(h/1000, LW_vs_h(:,i), 'LineWidth', 1.8, 'Color', C(i,:));
end
xlabel('Altitude, h (km)');
ylabel('L/W at \alpha^*');
title('L/W vs altitude (using \alpha^* that maximizes L/D over all \alpha and \DeltaT)');
legend(legtxt, 'Location','best');

% --- 2) T/W vs altitude for different Mach (use CT at same alpha*, dT* above) ---
TW_vs_h = zeros(numel(h), numel(Ma_list));
for i = 1:numel(Ma_list)
    TW_vs_h(:,i) = (q(:,i) .* S_ref .* CT_at_star_g(i)) ./ W;
end

figure('Name','T/W vs altitude (different Mach)'); hold on; grid on;
for i = 1:numel(Ma_list)
    plot(h/1000, TW_vs_h(:,i), 'LineWidth', 1.8, 'Color', C(i,:));
end
xlabel('Altitude, h (km)');
ylabel('T/W');
title('T/W vs altitude (for each Mach, using CT at \alpha^*, \DeltaT^* that gives global max L/D)');
legend(legtxt, 'Location','best');

% --- (Optional) show alpha* and dT* used for each Mach ---
disp('--- alpha* and dT* used (global max L/D over all alpha & dT) ---');
for i = 1:numel(Ma_list)
    fprintf('Ma=%.1f: max L/D=%.4f at alpha*=%.2f deg, dT*=%.2f, CL*=%.4f, CT*=%.4f\n', ...
        Ma_list(i), LD_max_global(i), alpha_star_g(i), dT_star_g(i), CL_at_star_g(i), CT_at_star_g(i));
end


%% ======================= Local function: ISA atmosphere =======================
function [T, a, p, rho] = isa_atm_simple(h)
% Simple ISA-like atmosphere (1976) to ~86 km.
% Input:
%   h   altitude [m] (vector)
% Output:
%   T   temperature [K]
%   a   speed of sound [m/s]
%   p   pressure [Pa]
%   rho density [kg/m^3]

h = max(h, 0);
R  = 287.05287;      % J/(kg*K)
g0 = 9.80665;        % m/s^2
gamma = 1.4;

% Define layers: base height hb [m], base temp Tb [K], base pressure pb [Pa], lapse Lb [K/m]
% Layers: 0-11km, 11-20km, 20-32km, 32-47km, 47-51km, 51-71km, 71-86km
hb = [0,    11000, 20000, 32000, 47000, 51000, 71000]';
Lb = [-0.0065, 0.0, 0.0010, 0.0028, 0.0, -0.0028, -0.0020]';

Tb = zeros(size(hb));
pb = zeros(size(hb));
Tb(1) = 288.15;         % sea level
pb(1) = 101325.0;       % sea level

% Precompute base T,p for each layer
for i = 2:numel(hb)
    h0 = hb(i-1);  T0 = Tb(i-1);  p0 = pb(i-1);  L0 = Lb(i-1);
    if abs(L0) < 1e-12
        Tb(i) = T0;
        pb(i) = p0 * exp(-g0*(hb(i)-h0)/(R*T0));
    else
        Tb(i) = T0 + L0*(hb(i)-h0);
        pb(i) = p0 * (Tb(i)/T0)^(-g0/(R*L0));
    end
end

T = zeros(size(h));
p = zeros(size(h));

for n = 1:numel(h)
    hn = h(n);

    % Find layer index
    idx = find(hn >= hb, 1, 'last');
    if idx == numel(hb) && hn > 86000
        idx = numel(hb); % clamp above top layer
        hn  = 86000;
    end

    h0 = hb(idx);  T0 = Tb(idx);  p0 = pb(idx);  L0 = Lb(idx);

    if abs(L0) < 1e-12
        T(n) = T0;
        p(n) = p0 * exp(-g0*(hn-h0)/(R*T0));
    else
        T(n) = T0 + L0*(hn-h0);
        p(n) = p0 * (T(n)/T0)^(-g0/(R*L0));
    end
end

rho = p ./ (R .* T);
a   = sqrt(gamma .* R .* T);
end