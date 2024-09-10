% function output = fun_rotm2quat(R)
% % converts rotation matrix to quaternion
% % output quaternion has nonnegative scalar component
% 
% tr = trace(R);
% 
% if tr > -1 && tr <= 3
%     q0 = (sqrt(1 + tr))/2;
%     vector_part = (fun_vex(R - R'))/(2*(sqrt(1 + tr)));
%     output = [q0; vector_part];
% 
% elseif tr == -1
%     [V, D] = eig(R);
%     % real eigenvalue at D(3,3)
%     v = V(:,3);
% 
%     output = [0; v/norm(v)];
% 
% else
%     error('Invalid rotation matrix')
% end


function q = fun_rotm2quat(R)
    % This function converts a rotation matrix to a quaternion.
    % Input: 
    %   R - 3x3 rotation matrix
    % Output:
    %   q - 4x1 quaternion [qw, qx, qy, qz]

    % Ensure the matrix is valid
    if ~isequal(size(R), [3, 3])
        error('Input matrix must be a 3x3 matrix');
    end
    
    % Pre-allocate the quaternion
    q = zeros(4, 1);

    % Compute the trace of the matrix
    tr = trace(R);

    % Depending on the trace, compute the quaternion
    if tr > 0
        S = sqrt(tr + 1.0) * 2; % S = 4 * qw
        q(1) = 0.25 * S;
        q(2) = (R(3, 2) - R(2, 3)) / S;
        q(3) = (R(1, 3) - R(3, 1)) / S;
        q(4) = (R(2, 1) - R(1, 2)) / S;
    elseif (R(1, 1) > R(2, 2)) && (R(1, 1) > R(3, 3))
        S = sqrt(1.0 + R(1, 1) - R(2, 2) - R(3, 3)) * 2; % S = 4 * qx
        q(1) = (R(3, 2) - R(2, 3)) / S;
        q(2) = 0.25 * S;
        q(3) = (R(1, 2) + R(2, 1)) / S;
        q(4) = (R(1, 3) + R(3, 1)) / S;
    elseif R(2, 2) > R(3, 3)
        S = sqrt(1.0 + R(2, 2) - R(1, 1) - R(3, 3)) * 2; % S = 4 * qy
        q(1) = (R(1, 3) - R(3, 1)) / S;
        q(2) = (R(1, 2) + R(2, 1)) / S;
        q(3) = 0.25 * S;
        q(4) = (R(2, 3) + R(3, 2)) / S;
    else
        S = sqrt(1.0 + R(3, 3) - R(1, 1) - R(2, 2)) * 2; % S = 4 * qz
        q(1) = (R(2, 1) - R(1, 2)) / S;
        q(2) = (R(1, 3) + R(3, 1)) / S;
        q(3) = (R(2, 3) + R(3, 2)) / S;
        q(4) = 0.25 * S;
    end

    % Normalize the quaternion
    q = q / norm(q);
end
