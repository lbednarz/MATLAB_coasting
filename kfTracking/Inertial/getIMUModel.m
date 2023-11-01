function [A, B, G] = getIMUModel(fb1, fb2, fb3, lat, lon, phi, theta, psi, ...
    phi_dot, theta_dot, psi_dot, settings)
% get jacobians before tracking starts
%E = [phi; theta; psi];
E_dot = [phi_dot; theta_dot; psi_dot]; 
Qbe_inv_sym = compute_Qbe_inv(phi, theta, psi);
R_BN_sym = enuToBody(phi, theta, psi);
w_IE_N = ecefToEnu(lat, lon)*[0;0;7.2722e-5];
f_b_star = [fb1;fb2;fb3];

dQbe_inv_dE = settings.dQ_be(phi, theta, psi);
dR_BN_dE = settings.dR_BN(phi, theta, psi);

% Computing coefficients of attitude vector 
w_NB_B = Qbe_inv_sym \ E_dot;  % As per inverse(Q_{be,inverse}) * dot(E)
w_IE_B = R_BN_sym * w_IE_N; % assuming w_EN_B = [0; 0; 0];
w_IB_B = w_IE_B + w_NB_B;
K = kron(eye(3), (w_IB_B - R_BN_sym*w_IE_B).') *  ...
        dQbe_inv_dE - kron(Qbe_inv_sym, w_IE_N.') * dR_BN_dE;
dE_N_coeff = [zeros(size(K)), K, zeros(size(K))];

% computing user velocity 
skewSymmetricMatrix = @(v) [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
f_B_N = R_BN_sym*f_b_star;
f_B_N_cross = skewSymmetricMatrix(f_B_N);
w_IE_N_cross = skewSymmetricMatrix(w_IE_N);
dv_N = [zeros(size(w_IE_N_cross)), w_IE_N_cross, f_B_N_cross];

% computing user position 
dr_N = eye(3);
dr_N_coeff = [zeros(size(dr_N)), dr_N, zeros(size(dr_N))];

% state transition matrix and control input matrix
A = [dr_N_coeff; dv_N; dE_N_coeff];
B = [zeros(size(dr_N)), zeros(size(dr_N));
     R_BN_sym, zeros(size(R_BN_sym));
     zeros(size(Qbe_inv_sym)), Qbe_inv_sym];
G = eye(9); % TODO 
end
