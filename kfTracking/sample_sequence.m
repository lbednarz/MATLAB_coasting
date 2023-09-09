% function to created sampled PRNs 
function sampled_sequence = sample_sequence(prn_sequence, sampling_rate)
    % Determine the number of samples per chip
    chip_rate = 1.023e6; % chips/1s 
    samples_per_chip = ceil(sampling_rate / chip_rate);
    samples_per_code = sampling_rate*1e-3; % 1 ms of data 
   %remainder = sampling_rate*1e-3 - samples_per_code; 

    % Repeat each element of the PRN sequence 'samples_per_chip' times
    sampled_sequence = repelem(prn_sequence, samples_per_chip);
    sampled_sequence = sampled_sequence(1:samples_per_code);
    
%     % encode into -1s and 1s
%     sampled_sequence(sampled_sequence == 1) = -1;
%     sampled_sequence(sampled_sequence == 0) = 1;
end