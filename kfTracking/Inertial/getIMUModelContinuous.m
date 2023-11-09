function [A, B, G] = getIMUModelContinuous(tau_acc, tau_gy)
% get jacobians before tracking starts
syms phi theta psi_var phi_dot theta_dot psi_var_dot fb1 fb2 fb3 lat lon;
E = [phi, theta, psi_var];
% Qbe_inv_sym = compute_Qbe_inv(phi, theta, psi_var);
R_BN_sym = enuToBodySA(phi, theta, psi_var);

%tensor = sym(zeros(3,3,3));
tensor2 = sym(zeros(3,3,3));
for i = 1:3
    for j = 1:3
        for k = 1:3
            %tensor(i,j,k) = diff(Qbe_inv_sym(i,j), E(k));
            tensor2(i,j,k) = diff(R_BN_sym(i,j), E(k));
        end
    end
end

% by-hand derivatives 
% dQbe_inv_dE = [0, cos(phi)*tan(theta), -sin(phi)*tan(theta);
%                0, sin(phi)*sec(theta)^2, cos(theta)*sec(theta)^2;
%                0, 0, 0;
%                0, -sin(phi), 0;
%                0, 0, 0;
%                0, 0, -cos(psi_var);
%                0, sin(psi_var)*sin(theta)/cos(theta)^2, cos(phi)*sin(theta)/cos(theta)^2;
%                0, cos(psi_var)/cos(theta), 0].';

%sm angle version
dQbe_inv_dE = [0, 0, 0;
               theta, phi, 0;
               -phi*theta, 1, 0;
               0, 0, 0;
               -phi, 0, 0;
               -1, 0, 0;
               0, 0, 0;
               1, phi*theta*theta, 0;
               -phi*theta, theta, 0];

%dQbe_inv_dE = reshape(tensor, 9, 3);
dR_BN_dE = reshape(tensor2, 9, 3);

E_dot       = [phi_dot; theta_dot; psi_var_dot]; 
Qbe_inv_sym = compute_Qbe_invSA(phi, theta, psi_var);
w_IE_N      = ecefToEnu(lat, lon)*[0;0;7.2722e-5];
f_b_star    = [fb1;fb2;fb3];

% Computing coefficients of attitude vector 
w_NB_B = Qbe_inv_sym \ E_dot;  % As per inverse(Q_{be,inverse}) * dot(E)
w_IE_B = R_BN_sym * w_IE_N; % assuming w_EN_B = [0; 0; 0];
w_IB_B = w_IE_B + w_NB_B;
K = kron(eye(3), (w_IB_B - R_BN_sym*w_IE_B).') *  ...
        dQbe_inv_dE - kron(Qbe_inv_sym, w_IE_N.') * dR_BN_dE;
dE_N_coeff = [zeros(3,3), zeros(3,3), K, ...
                zeros(3,3), -1*Qbe_inv_sym];

% computing user velocity 
skewSymmetricMatrix = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
f_B_N = R_BN_sym*f_b_star;
f_B_N_cross = skewSymmetricMatrix(f_B_N);
w_IE_N_cross = skewSymmetricMatrix(w_IE_N);
dv_N = [zeros(3,3), w_IE_N_cross, f_B_N_cross, ...
            -1*R_BN_sym, zeros(3,3)];

% computing user position 
dr_N = eye(3);
dr_N_coeff = [zeros(size(dr_N)), dr_N, zeros(size(dr_N)), ...
                zeros(size(dr_N)), zeros(size(dr_N))];

% accelerometer and gyro bias 
% skip 9 defined states and append to end
acc_coeff = -1\tau_acc;
gy_coeff = -1\tau_gy;
b_coeff = [zeros(6,9) diag([acc_coeff*ones(1,3), gy_coeff*ones(1,3)])];

% state transition matrix and control input matrix
% Set the tolerance value
tolerance = 1e-3;

A = vpa([dr_N_coeff; dv_N; dE_N_coeff; b_coeff],2);
B = vpa([zeros(size(dr_N)), zeros(size(dr_N));
     R_BN_sym, zeros(size(R_BN_sym));
     zeros(size(Qbe_inv_sym)), Qbe_inv_sym],2);
G = vpa([zeros(3,12); -1*R_BN_sym zeros(3,9); zeros(3,3) -1*Qbe_inv_sym ...
        zeros(3,6); zeros(6,6) eye(6)],2); 
end