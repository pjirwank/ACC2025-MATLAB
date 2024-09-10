%% Parameters

kp_bar = 1;

D = diag([1,2,3,4]); % eigenvalues of A

e1 = [1;0;0;0];
e2 = [0;1;0;0];
e3 = [0;0;1;0];
e4 = [0;0;0;1];

V = [e1 e2 e3 e4];

% check if all ei's are orthogonal
for i=1:1:3
    ei = V(:,i);
    for j=i+1:1:4
        ej = V(:,j);
        if abs(ei' * ej) > 0.0001
            error('Chosen vectors are not orthogonal.')
        end
    end
end

A = V * D * V';

% desired critical points
v1 = [1;0;0;0];
v2 = [1/sqrt(4);1/sqrt(4);1/sqrt(4);1/sqrt(4)];
v3 = [1/sqrt(4);-1/sqrt(4);1/sqrt(4);1/sqrt(4)];
v4 = [1/sqrt(4);1/sqrt(4);-1/sqrt(4);1/sqrt(4)];

% transformation matrix B such that B * vi = ei
B_inv = [v1 v2 v3 v4];
B = inv(B_inv);

%% Initializations

start_time = 0;
end_time = 10;
N = 100*(end_time - start_time) + 1;
time = linspace(start_time, end_time, N);
h = time(2) - time(1);

% initializations
R_init = eye(3);
Rhat_init = fun_axisangle(pi, [1;0;0]);

Ry_init = R_init;
Rtilde_init = Rhat_init' * R_init;
Rbar_init = Rhat_init * R_init';

% array initialiations
R_array = zeros(3,3,N);
R_array(:,:,1) = R_init;

Ry_array = zeros(3,3,N);
Ry_array(:,:,1) = Ry_init;

Rhat_array = zeros(3,3,N);
Rhat_array(:,:,1) = Rhat_init;

Rtilde_array = zeros(3,3,N);
Rtilde_array(:,:,1) = Rtilde_init;

Rbar_array = zeros(3,3,N);
Rbar_array(:,:,1) = Rbar_init;

Omega_array = zeros(3,N);
Omega_y_array = zeros(3,N);

for i=1:1:N
    Omega_array(:,i) = [sin(time(i)); cos(time(i)); sin(2*time(i))];
    Omega_y_array(:,i) = Omega_array(:,i);
end

% store morse function values for Rbar
morseFunc_array = zeros(N,1);
morseFunc_array(1) = fun_morseFuncOnRP3(fun_rotm2quat(Rbar_init), A, B);

potential_Rbar = zeros(N,1);
potential_Rbar(1) = fun_potential(Rbar_init);

%% main loop

for i=1:1:N-1
    % initialize
    R = R_array(:,:,i);
    Rhat = Rhat_array(:,:,i);
    Ry = Ry_array(:,:,i);
    Rtilde = Rtilde_array(:,:,i);
    Rbar = Rbar_array(:,:,i);
    Omega = Omega_array(:,i);
    Omega_y = Omega_y_array(:,i);
    
    % propagate
    R_next = fun_rotationPropagation(R, Omega, h);
    Rhat_next = fun_gradObserver(Rhat, Ry, Omega_y, A, B, kp_bar, h);
    
    % update
    R_array(:,:,i+1) = R_next;
    Rhat_array(:,:,i+1) = Rhat_next;
    Ry_array(:,:,i+1) = R_array(:,:,i+1);
    Rtilde_array(:,:,i+1) = Rhat_next' * R_next;
    Rbar_array(:,:,i+1) = Rhat_next * R_next';
    
    qbar_next = fun_rotm2quat(Rbar_array(:,:,i+1));
    morseFunc_array(i+1) = fun_morseFuncOnRP3(qbar_next, A, B);
    potential_Rbar(i+1) = fun_potential(Rbar_array(:,:,i+1));
end

%% plots
% 
figure(1)
plot(time, morseFunc_array - D(1,1)*ones(N, 1));
xlabel('Time $[s]$', 'Interpreter', 'latex')
ylabel('$f \circ \varphi (\bar{R}) - \lambda_1$', 'Interpreter', 'latex')
title('$f \circ \varphi (\bar{R}) - \lambda_1$', 'Interpreter', 'latex')
grid on

% figure(2)
% plot(time, potential_Rbar - fun_potential(fun_quat2rotm(v1)), 'LineWidth', 2)
% xlabel('Time $[s]$', 'Interpreter', 'latex')
% ylabel('$|\bar{R}_I|$', 'Interpreter', 'latex')
% title('$|\bar{R}_I|$', 'Interpreter', 'latex')
% grid on