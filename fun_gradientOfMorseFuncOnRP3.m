function output = fun_gradientOfMorseFuncOnRP3(v, A, B)
% v defines the point in RP3 at which the gradient is calculated
% A and B define the morse function on RP3
% output is the gradient

funcEval = fun_morseFuncOnRP3(v, A, B);

num = 2*((B' * A * B * v) - (funcEval * (B' * B) * v));
den = v' * (B' * B) * v;

output = num/den;
end