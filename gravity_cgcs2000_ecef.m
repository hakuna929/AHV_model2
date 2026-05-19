function g_ecef = gravity_cgcs2000_ecef(lat_deg, lon_deg, h)
%GRAVITY_CGCS2000_ECEF
% CGCS2000椭球下正常重力模型，输出ECEF重力加速度矢量 [m/s^2]

    % CGCS2000 / GRS80 椭球参数
    f  = 1 / 298.257222101;
    e2 = f * (2 - f);

    % Somigliana公式参数
    gamma_e = 9.7803267715;       % 赤道重力
    k = 0.00193185138639;

    lat = deg2rad(lat_deg);
    lon = deg2rad(lon_deg);

    sin_lat = sin(lat);

    % 椭球面正常重力
    gamma0 = gamma_e * (1 + k * sin_lat^2) / sqrt(1 - e2 * sin_lat^2);

    % 自由空气改正
    gamma_h = gamma0 - 3.086e-6 * h;

    % ECEF中"向上"方向
    up_ecef = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
    up_ecef = up_ecef / max(norm(up_ecef), 1e-12);

    % 重力矢量
    g_ecef = -gamma_h * up_ecef;
end