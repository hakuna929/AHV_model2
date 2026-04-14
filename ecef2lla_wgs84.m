function [lat, lon, h] = ecef2lla_wgs84(r_ecef)
    % WGS84 参数
    a = 6378137.0;
    f = 1 / 298.257223563;
    e2 = f * (2 - f);
    
    x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);
    
    % 经度可直接计算
    lon = atan2(y, x);
    
    % 纬度与高度迭代
    p = sqrt(x^2 + y^2);
    lat = atan2(z, p * (1 - e2));   % 初始近似
    for iter = 1:5
        sin_lat = sin(lat);
        N = a / sqrt(1 - e2 * sin_lat^2);
        h = p / cos(lat) - N;
        lat = atan2(z, p * (1 - e2 * N / (N + h)));
    end
end