% function [CL, CD, CY, Cl, Cm, Cn] = aero_coeffs(Ma, alpha, beta, de1, de2, dr, dT)
% % alpha, beta, de1, de2, dr 使用度要统一
% adeg = alpha*180/pi;
% bdeg = beta*180/pi;
% de1d = de1*180/pi;
% de2d = de2*180/pi;
% drd  = dr*180/pi;
% 
% 
% % C_La = C_La(Ma,α,δe1,δe2)
% CLa = 0.1498 - 0.02751*Ma + 0.07235*adeg ...
%      - 0.003368*adeg*Ma + 0.002343*Ma.^2 ...
%      + 0.001185*adeg.^2 + 0.006515*(de1d + de2d);
% 
% % C_Da = C_Da(Ma,α,δe1,δe2)
% CDa = 0.05099 - 0.004863*Ma + 0.002967*adeg ...
%      + 0.001364*adeg.^2 + 0.00053627*(de1d + de2d) ...
%      + 0.0001313*(de1d.^2 + de2d.^2);
% 
% % C_Y = C_Y(Ma,α,β,δr)  侧向力系数
% CY = 0.04626*bdeg - 0.002833*bdeg*Ma ...
%      - 0.0003691*adeg*bdeg ...
%      + (-0.007779 + 0.00050269*Ma)*drd;
% 
% %  C_l = C_l(Ma,β,δe1,δe2,δr) 滚转力矩系数
% Cl = (0.008739 - 0.00045027*Ma)*bdeg ...
%      + (0.003319 - 0.00021396*Ma)*drd ...
%      + (0.004906 - 0.00022831*Ma)*(de1d-de2d);
% 
% %  C_mα = C_mα(Ma,α,δe1,δe2)
% Cma = 0.2414 + 0.011618*Ma + 0.1012*adeg ...
%      + 0.00121*adeg.^2 ...
%      + (-0.03702 + 0.001733*Ma)*(de1d+de2d);
% 
% %  C_n = C_n(Ma,α,β,δr,δe1,δe2) 偏航力矩系数
% Cn = -0.001495*bdeg*Ma - 0.00028642*adeg*bdeg ...
%      + 0.01215*bdeg + (-0.01255 + 0.00098133*Ma)*drd ...
%      - 0.000027973*de1d.^3 - 0.0000787*(de1d.^2 - de2d.^2);
% 
% 
% 
% % C_Le = C_Le(Ma, alpha, beta, delta_T)
% CLe = 0.7215 ...
%     + 0.02635  * adeg ...
%     + 0.1147   * Ma*cosd(bdeg) ...
%     - 0.002795 * Ma*cosd(bdeg) * adeg ...
%     - 0.5782   * sqrt(Ma*cosd(bdeg)) ...
%     + 0.2894   * dT ...
%     - 0.004363 * adeg * dT ...
%     + 0.01083  * Ma*cosd(bdeg) * dT;
% 
% % C_De = C_De(Ma, alpha, beta, delta_T)
% CDe = 0.002339 * adeg ...
%     + 0.00012182 * adeg.^2 ...
%     - 0.00033126 * Ma*cosd(bdeg) * adeg ...
%     + 0.005557  * adeg * dT;
% 
% % C_me = C_me(Ma, alpha, beta, delta_T)
% Cme = -0.8297 ...
%     + 0.2703   * Ma*cosd(bdeg) ...
%     - 0.1133   * adeg ...
%     - 1.2315   * dT ...
%     + 0.01695  * adeg * dT ...
%     - 0.04602  * Ma*cosd(bdeg)* dT;
% 
% CL = CLa + CLe;
% CD = CDa + CDe;
% Cm = Cma + Cme;
% end


function [CL, CD, CY, Cl, Cm, Cn] = aero_coeffs(Ma, alpha, beta, de1, de2, dr, dT)
% 3DOF简化一致版气动系数
% 目的：
% 1) 3DOF不考虑舵偏动态影响 -> 舵偏输入不参与CL/CD/CY（固定为配平思想）
% 2) 去掉dT对气动系数耦合 -> 推进只由 thrust_coeffs 负责
% 3) 保留原接口，避免主程序改动
%
% 输入:
%   Ma, alpha(rad), beta(rad), de1,de2,dr,dT  (后四者在本模型中不参与气动力计算)
% 输出:
%   CL,CD,CY 参与3DOF平动；Cl,Cm,Cn仅占位返回0

    %#ok<*INUSD>  % de1,de2,dr,dT在此简化模型中不使用

    %% -------- 单位转换 --------
    adeg = alpha * 180/pi;
    bdeg = beta  * 180/pi;

    %% -------- 输入限幅（数值稳定）--------
    Ma   = min(max(Ma,   0.0), 8.0);
    adeg = min(max(adeg,-10.0), 10.0);
    bdeg = min(max(bdeg,-6.0), 6.0);

    %% -------- 升力系数 CL（仅Ma,alpha,beta）--------
    % 采用你原模型去舵偏后的主体项：
    CL = 0.1498 ...
       - 0.02751 * Ma ...
       + 0.07235 * adeg ...
       - 0.003368 * adeg * Ma ...
       + 0.002343 * Ma.^2 ...
       + 0.001185 * adeg.^2;

    % 轻微侧滑惩罚（可关）
    CL = CL - 0.0005 * bdeg.^2;

    %% -------- 阻力系数 CD（仅Ma,alpha,beta）--------
    CD = 0.05099 ...
       - 0.004863 * Ma ...
       + 0.002967 * adeg ...
       + 0.001364 * adeg.^2;

    % 侧滑附加阻力（可调小）
    CD = CD + 0.00025 * bdeg.^2;

    %% -------- 侧力系数 CY（简化）--------
    % 若你希望"纯纵向/平面"先收敛，可直接 CY=0
    CY = 0.04626 * bdeg ...
       - 0.002833 * bdeg * Ma ...
       - 0.0003691 * adeg * bdeg;
    % 如果要完全关闭侧向力，改成：
    % CY = 0;

    %% -------- 系数硬限幅（防数值炸裂）--------
    CL = min(max(CL, -2.2), 2.2);
    CD = min(max(CD,  0.015), 1.2);
    CY = min(max(CY, -0.8), 0.8);

    %% -------- 3DOF中不用的力矩系数 --------
    Cl = 0.0;
    Cm = 0.0;
    Cn = 0.0;
end