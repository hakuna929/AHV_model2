% 绘制密度随高度变化 平滑后的结果

clc;
clear;
close all;

%% 高度范围
h = linspace(0,50e3,2000);

rho_raw = zeros(size(h));
rho_smooth = zeros(size(h));

%% 计算
for i = 1:length(h)



    [~,rho_smooth(i)] = atmos_simple(h(i));

end

%% 绘图
figure;


plot(h/1000,rho_smooth,...
    'LineWidth',2);

grid on;

xlabel('Altitude (km)');
ylabel('Density (kg/m^3)');

title('Atmospheric Density');


xline(32.1619,'--');

set(gca,'FontSize',12);

%% 局部放大
figure;

idx = h/1000 > 28 & h/1000 < 36;


plot(h(idx)/1000,rho_smooth(idx),...
    'LineWidth',2);

grid on;

xlabel('Altitude (km)');
ylabel('Density (kg/m^3)');

title('Density Near 32 km Transition');


xline(32.1619,'--');

set(gca,'FontSize',12);
