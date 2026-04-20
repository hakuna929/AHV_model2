function r_ecef = lla2ecef_cgcs2000(lat, lon, alt)
% LLA2ECEF_CGCS2000 将CGCS2000椭球下的经纬高转换为ECEF直角坐标
%
%  输入：
%      lat : 纬度 [deg]，北纬为正 (-90 ~ 90)
%      lon : 经度 [deg]，东经为正 (-180 ~ 180)
%      alt : 海拔高度 [m]
%
%  输出：
%      r_ecef : 3x1 矢量，ECEF坐标系下的位置 [m] （原点为地心）
%
%  说明：
%      采用CGCS2000椭球参数 (a = 6378137 m, f = 1/298.257222101)

% CGCS2000 椭球定义常数
a  = 6378137.0;           % 长半轴 [m]
f  = 1 / 298.257222101;   % 扁率
e2 = f * (2 - f);         % 第一偏心率平方

% 角度转弧度
lat_rad = deg2rad(lat);
lon_rad = deg2rad(lon);

% 计算卯酉圈曲率半径 N
sin_lat = sin(lat_rad);
N = a ./ sqrt(1 - e2 * sin_lat.^2);

% 计算地心直角坐标分量
x = (N + alt) .* cos(lat_rad) .* cos(lon_rad);
y = (N + alt) .* cos(lat_rad) .* sin(lon_rad);
z = (N * (1 - e2) + alt) .* sin(lat_rad);

% 输出列矢量
r_ecef = [x; y; z];
end