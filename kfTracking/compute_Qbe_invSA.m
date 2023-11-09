function Qbe_inv = compute_Qbe_invSA(phi, theta, psi)
    % Computes the Q_{be,inverse} matrix based on the given Euler angles

    Qbe_inv = [1, phi*theta, theta;
               0, 1, -phi;
               0, phi, 1];
end