function output = fun_quatAxisAngle(theta, axis)

output = zeros(4,1);
output(1) = cos(theta/2);
output(2:4) = sin(theta/2)*axis;
end