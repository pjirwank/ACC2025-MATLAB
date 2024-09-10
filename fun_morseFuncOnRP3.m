function output = fun_morseFuncOnRP3(v, A, B)
% outputs f([v]) = <v, Av>/<v,v>
% A = V * D * inv(V)

output = (v' * B' * A * B * v)/(v' * (B' * B) * v);
end
