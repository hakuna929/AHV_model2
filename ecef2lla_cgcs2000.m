function [lat_deg, lon_deg, h] = ecef2lla_cgcs2000(r_ecef)
% ECEF2LLA_CGCS2000  将CGCS2000椭球下的ECEF直角坐标转换为经纬高(LLA)
% 输出：
%   lat_deg, lon_deg: [deg]
%   h: [m]

    % CGCS2000 参数（与GRS80一致）
    a = 6378137.0;
    f = 1 / 298.257222101;
    e2 = f * (2 - f);

    x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);

    % 经度（弧度）
    lon = atan2(y, x);

    % 迭代求纬度/高度（弧度迭代）
    p = hypot(x, y);
    lat = atan2(z, p * (1 - e2));  % 初值(弧度)

    for iter = 1:7
        sin_lat = sin(lat);
        N = a / sqrt(1 - e2 * sin_lat^2);
        h = p / max(cos(lat), 1e-12) - N;
        lat = atan2(z, p * (1 - e2 * N / (N + h)));
    end

    % 输出转成度
    lat_deg = rad2deg(lat);
    lon_deg = rad2deg(lon);
end