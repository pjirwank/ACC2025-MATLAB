function output = fun_potential(R)
% potential function
    output = sqrt((1/4)*trace(eye(3) - R));
end