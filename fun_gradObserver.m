function output = fun_gradObserver(Rhat, Ry, Omega_y, A, B, kp_bar, h)
% output = Rhat_plus % estimate at next time instant

Rbar_y = Rhat * Ry';

qbar_y = fun_rotm2quat(Rbar_y); % quaternion corresponding to Rbar
grad = fun_gradientOfMorseFuncOnRP3(qbar_y, A, B);
pullback_grad = fun_inverseTangentMapVarphi(qbar_y, grad);

Omega_Rhat = Omega_y - kp_bar * fun_vex(Rhat' * pullback_grad * Ry);
output = fun_rotationPropagation(Rhat, Omega_Rhat, h);
end