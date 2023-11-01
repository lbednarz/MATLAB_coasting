function Qbe_inv = compute_Qbe_inv(phi, theta, psi)
    % Computes the Q_{be,inverse} matrix based on the given Euler angles

    Qbe_inv = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
               0, cos(phi), -sin(phi);
               0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
end
