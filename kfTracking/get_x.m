% function that accepts a time of travel, PRN, and sample rate it returns 
% the PRN sequence consistent with that state 
function x_tau_seq = get_x(rem_tot, PRN, Ts)
    % Generate CA code
    caCode = generateCAcode(PRN);

    % Sample the sequence
    sampled_sequence = sample_sequence(caCode, Ts);

    % Calculate the shift index
    Tc = 1e-3; % code period is 1 ms
    index = floor(rem_tot / Tc * length(sampled_sequence));

    % Shift the sequence
    x_tau_seq = circshift(sampled_sequence, [-index, 0]);
end