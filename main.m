clear;
clc;
close all;

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
maxNoisyAngvelMagnitude = (10)*(pi/180); % rad/s

%% Initializations
start_time = 0;
end_time = 15;
N = 500*(end_time - start_time) + 1;
time = linspace(start_time, end_time, N);
h = time(2) - time(1);


axis_EKF = [0;1;0];
Rhat_init = fun_axisangle(180.*pi/180, axis_EKF);
R_init = eye(3);
Rbar_init = Rhat_init * R_init';
Rtilde_init = Rhat_init' * R_init;
q_init      = 0;

% array initializations
q_array = zeros(N,1);
q_array(1) = q_init;

R_array = zeros(3,3,N);
R_array(:,:,1) = R_init;

Rhat_array = zeros(3,3,N); 
Rhat_array(:,:,1) = Rhat_init;

Rhat_PCF_array = zeros(3,3,N);
Rhat_PCF_array(:,:,1) = Rhat_init;

Rhat_EKF_array = zeros(3,3,N);
Rhat_EKF_array(:,:,1) = Rhat_init;

Rtilde_array = zeros(3,3,N);
Rtilde_array(:,:,1) = Rtilde_init;

Rtilde_PCF_array = zeros(3,3,N);
Rtilde_PCF_array(:,:,1) = Rtilde_init;

Rtilde_EKF_array = zeros(3,3,N);
Rtilde_EKF_array(:,:,1) = Rtilde_init;

Rbar_array = zeros(3,3,N);
Rbar_array(:,:,1) = Rbar_init;

Rbar_PCF_array = zeros(3,3,N);
Rbar_PCF_array(:,:,1) = Rbar_init;

Rbar_EKF_array = zeros(3,3,N);
Rbar_EKF_array(:,:,1) = Rbar_init;

Omega_array = zeros(3,N); 
for i=1:1:N
    Omega_array(:,i) = [sin(time(i)), cos(time(i)), sin(2*time(i))];
end

jumps = zeros(N,1); % number of jumps
jump_times = [];
jump_index = [];

potential_Rbar = zeros(N,1);
potential_Rbar(1) = fun_potential(Rbar_init);

potential_Rbar_PCF = zeros(N,1);
potential_Rbar_PCF(1) = fun_potential(Rbar_init);

potential_Rbar_EKF = zeros(N,1);
potential_Rbar_EKF(1) = fun_potential(Rbar_init);

P_array = zeros(3,3,N);
P_array(:,:,1) = 1e-5 * eye(3);

cov_Ry = 1e-2 * eye(3);

for i = 1:1:N-1

    % if noise == 0
    %     noise_angvel_EKF = zeros(3,1);
    %     noise_rotm_EKF = eye(3);
    % elseif noise == 1
    %     direction_angvel_noise = rand(3,1);
    %     if norm(direction_angvel_noise) < 0.001
    %         direction_angvel_noise = rand(3,1);
    %     end
    %     direction_angvel_noise = normalize(direction_angvel_noise);
    %     max_noise_angvel = 10*pi/180; % rad/s
    %     noise_angvel_EKF = (-max_noise_angvel + 2*max_noise_angvel*rand())*direction_angvel_noise;
    %     max_noise_angle = 5*pi/180; %rad
    %     noise_angle_EKF = (-max_noise_angle + 2*max_noise_angle*rand());
    %     axis_noise_rotm = rand(3,1);
    %     if norm(axis_noise_rotm) < 0.001
    %         axis_noise_rotm = rand(3,1);
    %     end
    %     axis_noise_rotm = axis_noise_rotm/norm(axis_noise_rotm);
    %     noise_rotm_EKF = fun_axisangle(noise_angle_EKF, axis_noise_rotm);
    % else
    %     error('Incorrect noise value; not 0 or 1.')
    % end
    
    if noise == 0
        noise_angvel_PCF = zeros(3,1);
        noise_rotm_PCF = eye(3);
        noise_angvel_EKF = zeros(3,1);
        noise_rotm_EKF = eye(3);
    elseif noise == 1
        % noise for PCF and hybrid
        Rtilde_PCF_axis = fun_findaxis(Rtilde_PCF_array(:,:,i));
        Rbar_PCF_angle = real((acos(1 - 2*(potential_Rbar_PCF(i))))); %*(180/pi);
        
        if fun_potential(Rtilde_PCF_array(:,:,i)) > 0.9999
%             fprintf('case 1 \n')
            axis_Rtilde_PCF = fun_findaxis(Rtilde_PCF_array(:,:,i));
            sinTheta = sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i))));
            noise_angvel_PCF = maxNoisyAngvelMagnitude*sinTheta * axis_Rtilde_PCF;
        else
%             fprintf('case 2\n')
            axis_Rtilde_PCF = fun_findaxis(Rtilde_PCF_array(:,:,i));
            sinTheta = sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i))));
            noise_angvel_PCF = -maxNoisyAngvelMagnitude*sign(sinTheta) * axis_Rtilde_PCF;
        end

        % if potential_Rbar_EKF(i) < 0.9999
        %     noise_angvel_PCF = -sign(sin(Rbar_PCF_angle))*Rtilde_PCF_axis;
        %     % noise_angvel_PCF = -sign(sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i)))))*Rtilde_PCF_axis;
        % else
        %     % noise_angvel_PCF = sin(Rbar_PCF_angle)*Rtilde_PCF_axis;
        %     noise_angvel_PCF = (sin(2*asin(fun_potential(Rtilde_PCF_array(:,:,i)))))*Rtilde_PCF_axis;
        % 
        % end

        max_noise_angle_PCF = 10*pi/180; %rad
        noise_angle_PCF = (-max_noise_angle_PCF + 2*max_noise_angle_PCF*rand());
        axis_noise_rotm_PCF = rand(3,1);
        if norm(axis_noise_rotm_PCF) < 0.001
            axis_noise_rotm_PCF = rand(3,1);
        end
        axis_noise_rotm_PCF = axis_noise_rotm_PCF/norm(axis_noise_rotm_PCF);
        noise_rotm_PCF = fun_axisangle(noise_angle_PCF, axis_noise_rotm_PCF);
        
        
        % noise for EKF
        Rtilde_EKF_axis = fun_findaxis(Rtilde_EKF_array(:,:,i));
        Rbar_EKF_angle = real((acos(1 - 2*(potential_Rbar_EKF(i))))); %*(180/pi);
        if potential_Rbar_EKF(i) < 0.9999
            noise_angvel_EKF = -sign(sin(Rbar_EKF_angle))*Rtilde_EKF_axis;
        else
            noise_angvel_EKF = sin(Rbar_EKF_angle)*Rtilde_EKF_axis;
        end

        max_noise_angle = 10*pi/180; %rad
        noise_angle_EKF = (-max_noise_angle + 2*max_noise_angle*rand());
        axis_noise_rotm = rand(3,1);
        if norm(axis_noise_rotm) < 0.001
            axis_noise_rotm = rand(3,1);
        end
        axis_noise_rotm = axis_noise_rotm/norm(axis_noise_rotm);
        noise_rotm_EKF = fun_axisangle(noise_angle_EKF, axis_noise_rotm);
    else
        error('Incorrect noise value; not 0 or 1.')
    end
    
    Ry_PCF = noise_rotm_PCF * R_array(:,:,i);
    Ry_EKF =  noise_rotm_EKF * R_array(:,:,i);
    Omega_y_PCF = Omega_array(:,i) + noise_angvel_PCF;
    Omega_y_EKF = Omega_array(:,i) + noise_angvel_EKF;

    R_array(:,:,i+1) = fun_rotationPropagation(R_array(:,:,i), Omega_array(:,i), h);
    [Rhat_array(:,:,i+1), q_array(i+1), jumps(i+1)] = fun_hybridPCF(q_array(i), Rhat_array(:,:,i), Ry_PCF, Omega_y_PCF, A, B, kp, kp_bar, c0, c1, h, jumps(i));
    Rhat_PCF_array(:,:,i+1) = fun_passiveComplementaryFilter(Rhat_PCF_array(:,:,i), Ry_PCF, Omega_y_PCF, kp, h);
    [Rhat_EKF_array(:,:,i+1), P_array(:,:,i+1)] = fun_EKF(Rhat_EKF_array(:,:,i), P_array(:,:,i), Ry_EKF, Omega_y_EKF, cov_Ry, h);

    Rbar_array(:,:,i+1) = Rhat_array(:,:,i+1) * R_array(:,:,i+1)';
    Rtilde_array(:,:,i+1) = Rhat_array(:,:,i+1)' * R_array(:,:,i+1);
    potential_Rbar(i+1) = fun_potential(Rbar_array(:,:,i+1));

    Rbar_PCF_array(:,:,i+1) = Rhat_PCF_array(:,:,i+1) * R_array(:,:,i+1)';
    Rtilde_PCF_array(:,:,i+1) = Rhat_PCF_array(:,:,i+1)' * R_array(:,:,i+1);
    potential_Rbar_PCF(i+1) = fun_potential(Rbar_PCF_array(:,:,i+1));
    
    Rbar_EKF_array(:,:,i+1) = Rhat_EKF_array(:,:,i+1) * R_array(:,:,i+1)';
    Rtilde_EKF_array(:,:,i+1) = Rhat_EKF_array(:,:,i+1)' * R_array(:,:,i+1);
    potential_Rbar_EKF(i+1) = fun_potential(Rbar_EKF_array(:,:,i+1));

    % calculate jump times
    if jumps(i+1) - jumps(i) ~= 0
        jump_times(end+1) = time(i);
        jump_index(end+1) = i;
    end
end

number_of_jumps = length(jump_index);

%% plots

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
hold on
plot(time, potential_Rbar_EKF, '--m', 'LineWidth', 2)
ax = gca;
ax.FontSize = 20;
legend('Hybrid Filter', '', '', '', '', '', '', 'PCF', 'EKF', 'Interpreter', 'latex', 'FontSize', 35)
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