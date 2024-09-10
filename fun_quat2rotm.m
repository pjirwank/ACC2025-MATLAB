function output = fun_quat2rotm(q)
% outputs rotation matrix correspondng to quaternion q

n = q(1);
e = q(2:4);

output = eye(3) + 2*n*fun_cross(e) + 2*(fun_cross(e))^2;

end