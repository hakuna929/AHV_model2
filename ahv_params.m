function p = ahv_params()
% 飞行器与环境参数

p.m   = 671.33;      % 结构质量 (kg)，不含燃油
p.mf  = 100;         % 初始燃油 (kg)

p.L   = 4.27;        % 机体长度 (m)
p.S   = 0.2986;      % 参考面积 (m^2)
p.c   = 0.3732;      % 平均气动弦长 c_bar (m)
p.b   = 0.8;         % 翼展参考 (m)

p.Ix  = 34.13;       % kg m^2
p.Iy  = 1040;
p.Iz  = 1034;
p.Ixz = 430;

% 一些参考量，用于推进系数 -> 实际推力
p.H_ref  = 30e3;     % m
p.Ma_ref = 6;        % 参考马赫
[a_ref, rho_ref] = atmos_simple(p.H_ref);
p.rho_ref = rho_ref;
p.V_ref   = p.Ma_ref * a_ref;
end