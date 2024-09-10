%% Parameters

theta_c0 = (120)*(pi/180);
theta_c1 = (60)*(pi/180);
c0 = fun_potential(fun_axisangle(theta_c0, [1;0;0]));
c1 = fun_potential(fun_axisangle(theta_c1, [1;0;0]));

if (c1 > c0)
    error('Invalid values for c0 and c1.')
end

kp = 1;
kp_bar = 1;

D = diag([1,2,3,4]); % eigenvalues of A

% eigenvectors on A
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
            error('Chosen vectors are not orthogonal for i = %i and j = %i.', i, j);
        end
    end
end

A = V * D * V';

% desired critical points
% set such that \varphi^{-1}([v_i]) has rotation angle less than theta_c1
theta1 = (30)*(pi/180);
theta2 = (45)*(pi/180);
theta3 = (45)*(pi/180);
theta4 = (45)*(pi/180);

axis1 = [1;0;0];
axis2 = [0;1;0];
axis3 = [0;0;1];
axis4 = sqrt(1/3)*[1;1;1];

v1 = fun_quatAxisAngle(theta1, axis1);
v2 = fun_quatAxisAngle(theta2, axis2);
v3 = fun_quatAxisAngle(theta3, axis3);
v4 = fun_quatAxisAngle(theta4, axis4);

% v1 = [1;0;0;0];
% v2 = [1/sqrt(4);1/sqrt(4);1/sqrt(4);1/sqrt(4)];
% v3 = [1/sqrt(4);-1/sqrt(4);1/sqrt(4);1/sqrt(4)];
% v4 = [1/sqrt(4);1/sqrt(4);-1/sqrt(4);1/sqrt(4)];

mat = [v1 v2 v3 v4]; % matrix of the above vectors

for i=1:1:4
    vi = mat(:,i);
    vi = vi/norm(vi);
    if fun_potential(fun_quat2rotm(vi)) > c1
        error("The critical point i = %i not in the interior of the set D1.", i)
    end
end

% transformation matrix B such that B * vi = ei
B_inv = [v1 v2 v3 v4];
B = inv(B_inv);

%% Noise parameters

noise = 1; % 0 = no noise, 1 = noisy measurements

maxNoisyRotationAngle = (10)*(pi/180); % rad
maxNoisyAngvelMagnitude = 0.05; % rad/s

noisyAxis_init = rand(3,1);
noisyAxis_init = noisyAxis_init/norm(noisyAxis_init); % unit vector
noisyAngle_init = -maxNoisyRotationAngle + 2*maxNoisyRotationAngle*rand();

if noise == 0
    noisyRy_init = eye(3);
    noisyAngvel_init = zeros(3,1);
elseif noise == 1
    noisyRy_init = fun_axisangle(noisyAngle_init, noisyAxis_init);
    noisyAngvel_direction = rand(3,1);
    noisyAngvel_direction = noisyAngvel_direction/norm(noisyAngvel_direction);
    noisyAngvel_magnitude = -maxNoisyAngvelMagnitude + 2*maxNoisyAngvelMagnitude*rand();
    noisyAngvel_init = noisyAngvel_magnitude * noisyAngvel_direction;
else
    error('Wrong input for "noise".')
end

%% Initializations

start_time = 0;
end_time = 20;
N = 500*(end_time - start_time) + 1;
time = linspace(start_time, end_time, N);
h = time(2) - time(1);

% initializations
R_init      = eye(3);
Rhat_init   = fun_axisangle((180)*(pi/180), [1;0;0]);
Rtilde_init = Rhat_init' * R_init;
Rbar_init   = Rhat_init * R_init';
q_init      = 0;

% array initialiations
q_array = zeros(N,1);
q_array(1) = q_init;

R_array = zeros(3,3,N);
R_array(:,:,1) = R_init;

Ry_array = zeros(3,3,N);
% Ry_array(:,:,1) = Ry_init;

Rhat_array = zeros(3,3,N); 
Rhat_array(:,:,1) = Rhat_init;

Rhat_PCF_array = zeros(3,3,N);
Rhat_PCF_array(:,:,1) = Rhat_init;

Rtilde_array = zeros(3,3,N);
Rtilde_array(:,:,1) = Rtilde_init;

Rtilde_PCF_array = zeros(3,3,N);
Rtilde_PCF_array(:,:,1) = Rtilde_init;

Rbar_array = zeros(3,3,N);
Rbar_array(:,:,1) = Rbar_init;

Rbar_PCF_array = zeros(3,3,N);
Rbar_PCF_array(:,:,1) = Rbar_init;

Omega_array = zeros(3,N); % true ang vel
for i=1:1:N
    Omega_array(:,i) = [sin(time(i)); cos(time(i)); sin(2*time(i))];
end

Omega_y_array = zeros(3,N); % measured ang vel
% Omega_y_array(:,1) = Omega_array(:,1) + noisyAngvel_init;

jumps = zeros(N,1); % number of jumps
jump_times = [];
jump_index = [];

potential_Rbar = zeros(N,1);
potential_Rbar(1) = fun_potential(Rbar_init);

potential_Rbar_PCF = zeros(N,1);
potential_Rbar_PCF(1) = fun_potential(Rbar_init);

%% Main loop

for i=1:1:N-1
    
    % Noise calculation
    if noise == 0
        noisyRy = eye(3);
        noisyAngvel = zeros(3,1);
    elseif noise == 1
        
        % set attitude measurement noise
        noisyAxis = rand(3,1);
        noisyAxis = noisyAxis/norm(noisyAxis); % unit vector
        noisyAngle = -maxNoisyRotationAngle + 2*maxNoisyRotationAngle*rand();
        noisyRy = fun_axisangle(noisyAngle, noisyAxis);
        
        % set gyroscope measurement noise
        if fun_potential(Rtilde_PCF_array(:,:,i)) > 0.9999
%             fprintf('case 1 \n')
            axis_Rtilde_PCF = fun_findaxis(Rtilde_PCF_array(:,:,i));
            sinTheta = sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i))));
            noisyAngvel = maxNoisyAngvelMagnitude*sinTheta * axis_Rtilde_PCF;
        else
%             fprintf('case 2\n')
            axis_Rtilde_PCF = fun_findaxis(Rtilde_PCF_array(:,:,i));
            sinTheta = sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i))));
            noisyAngvel = -maxNoisyAngvelMagnitude*sign(sinTheta) * axis_Rtilde_PCF;
        end
        
%         noisyAngvel_direction = rand(3,1);
%         noisyAngvel_direction = noisyAngvel_direction/norm(noisyAngvel_direction);
%         noisyAngvel_magnitude = -maxNoisyAngvelMagnitude + 2*maxNoisyAngvelMagnitude*rand();
%         noisyAngvel = noisyAngvel_magnitude * noisyAngvel_direction;
    else
        error('Wrong input for "noise".')
    end
    
    % Variable initializations
    q       = q_array(i);
    R       = R_array(:,:,i);
    Rhat    = Rhat_array(:,:,i);
    Ry      = noisyRy * R;
    Omega   = Omega_array(:,i);
    Omega_y = Omega + noisyAngvel;
    
    Rhat_PCF = Rhat_PCF_array(:,:,i);
    
    Ry_array(:,:,i) = Ry;
    Omega_y_array(:,i) = Omega_y;
    
    totalJumps = jumps(i);
        
    % Propagation
    R_next = fun_rotationPropagation(R, Omega, h);
    [Rhat_next, q_next, jump_next] = fun_hybridPCF(q, Rhat, Ry, Omega_y, A, B, kp, kp_bar, c0, c1, h, totalJumps);
    Rhat_PCF_next = fun_passiveComplementaryFilter(Rhat_PCF, Ry, Omega_y, kp, h);
    
    Rtilde_next = Rhat_next' * R_next;
    Rbar_next = Rhat_next * R_next';
    Rtilde_PCF_next = Rhat_PCF_next' * R_next;
    Rbar_PCF_next = Rhat_PCF_next * R_next';
    
    % Update array values
    q_array(i+1) = q_next;
    R_array(:,:,i+1) = R_next;
    Rhat_array(:,:,i+1) = Rhat_next;
    Rhat_PCF_array(:,:,i+1) = Rhat_PCF_next;
    Rtilde_array(:,:,i+1) = Rtilde_next;
    Rbar_array(:,:,i+1) = Rbar_next;
    Rbar_PCF_array(:,:,i+1) = Rbar_PCF_next;
    Rtilde_PCF_array(:,:,i+1) = Rtilde_PCF_next;
    jumps(i+1) = jump_next;
    
    potential_Rbar(i+1) = fun_potential(Rbar_next);
    potential_Rbar_PCF(i+1) = fun_potential(Rbar_PCF_next);
    
    % calculate jump times
    if jumps(i+1) - jumps(i) ~= 0
        jump_times(end+1) = time(i);
        jump_index(end+1) = i;
    end
end

Ry_array(:,:,end) = R_array(:,:,end);
Omega_y_array(:,end) = Omega_array(:,end);

number_of_jumps = length(jump_index);

%% Plots

clf;
close all;

figure(1)
plot(time, potential_Rbar, 'LineWidth', 2)
if number_of_jumps > 0
    for i=1:1:number_of_jumps
        ind = jump_index(i);
        hold on;
        plot(time(ind), potential_Rbar(ind), 'color', 'red', 'marker','x', 'linewidth', 2, 'MarkerSize',12);
        hold on;
        plot(time(ind+1), potential_Rbar(ind+1), 'color', 'red', 'marker','o', 'linewidth', 1, 'MarkerSize',9);
        hold on;
        plot([time(ind) time(ind+1)], [potential_Rbar(ind) potential_Rbar(ind+1)], '--r', 'LineWidth', 2);
    end
end
hold on
plot(time, potential_Rbar_PCF, 'color', 'black', 'LineWidth', 2)
ax = gca;
ax.FontSize = 20;
legend('Hybrid Filter', '', '', '', '', '', '', 'PCF', 'Interpreter', 'latex', 'FontSize', 35)
xlabel("$t \: [s]$", 'Interpreter', 'latex', 'FontSize', 35)
ylabel({'$|\bar{R}|_I $'}, 'interpreter', 'latex', 'FontSize', 35)
ylim([0, 1.1])
grid on

figure(2)
plot(time, q_array, 'o')
if number_of_jumps == 0
    return;
else
    for i=1:1:number_of_jumps
        ind = jump_index(i);
        hold on;
        plot(time(ind), q_array(ind), 'color', 'red', 'marker','x', 'linewidth', 2, 'MarkerSize',12);
        hold on;
        plot(time(ind+1), q_array(ind+1), 'color', 'red', 'marker','o', 'linewidth', 1, 'MarkerSize',9);
        hold on;
        plot([time(ind) time(ind+1)], [q_array(ind) q_array(ind+1)], '--r', 'LineWidth', 2);
    end
end
ax = gca;
ax.FontSize = 20;
xlabel("$t \: [s]$", 'Interpreter', 'latex', 'FontSize', 35)
ylabel({'$q$'}, 'interpreter', 'latex', 'FontSize', 35)
yticks([0,1])
grid on