function [value, isterminal, direction] = event_reach_target(t, X, rt_n, R_hit)
% 当与目标距离小于 R_hit 时终止积分
rm_n = X(1:3);
R = rt_n - rm_n;
range = norm(R);

value = range - R_hit;   % value=0 时触发
isterminal = 1;          % 终止
direction = -1;          % 只在从大到小穿越时触发
end