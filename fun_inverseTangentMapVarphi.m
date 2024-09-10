function output = fun_inverseTangentMapVarphi(q, v)
% v\in T_[q] RP3
% output \in T_R SO(3)

n = q(1);
e = q(2:4);
R = fun_quat2rotm(q);

ang_vel = 2 * [-e'; n*eye(3)+fun_cross(e)]' * v;

output = R*fun_cross(ang_vel);
end
