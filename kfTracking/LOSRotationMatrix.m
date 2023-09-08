function R = LOSRotationMatrix(r_receiver, r_SV)
    % Compute LOS vector
    zb = r_SV - r_receiver;
    zb = zb / norm(zb);
    
    % Earth's rotation axis in ECEF
    e_rotation = [0; 0; 1];
    
    % Compute yb
    yb = cross(zb, e_rotation);
    yb = yb / norm(yb);
    
    % Compute xb
    xb = cross(yb, zb);
    
    % Form the rotation matrix
    R = [xb, yb, zb];
end
