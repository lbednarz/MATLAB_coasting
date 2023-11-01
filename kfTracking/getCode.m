function x_tau_seq = getCode(rem_tot, PRN, Ts, fd)
    % Generate CA code
    caCode = generateCAcode(PRN);

    % Calculate the effective code frequency
    fc = 1.023e6; % For GPS C/A code
    feff = fc * (1 + fd/fc); 

    % Calculate resampling factor
    resample_factor = fc / feff;

    % Resample the sequence
    % factor of 100 comes from "resample" requiring integer aguments
    resampled_sequence = resample(caCode, floor(resample_factor * 100), 100);

    % Sample the resampled sequence
    sampled_sequence = sample_sequence(resampled_sequence, Ts);

    % Calculate the shift index
    Tc = 1; % code period is 1 ms
    index = floor(rem_tot / Tc * length(sampled_sequence));

    % Shift the sequence
    x_tau_seq = circshift(sampled_sequence, -index);
end
