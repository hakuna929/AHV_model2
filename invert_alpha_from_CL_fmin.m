function alpha_deg = invert_alpha_from_CL_fmin(Ma, dT, CL_req, alpha_min, alpha_max)
    if nargin < 4, alpha_min = -2; end
    if nargin < 5, alpha_max = 10; end

    obj = @(a) (aero_coeffs(Ma, a, dT) - CL_req).^2;
    alpha_deg = fminbnd(obj, alpha_min, alpha_max);

    alpha_deg = min(max(alpha_deg, alpha_min), alpha_max);
end