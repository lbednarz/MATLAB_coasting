function R_NB = enuToBodySA(phi, theta, psi)
    R_NB = [1, theta*phi + psi, -1*theta + psi*phi;
            -psi, -psi*theta*phi + 1, psi*theta + phi;
            theta, -phi, 1];
end