% function that accepts a time of travel, PRN, and sample rate it returns 
% the PRN sequence consistent with that state 
function x_tau_seq = get_x(rem_tot, PRN, Ts)
    % Generate CA code
    caCode = generateCAcode(PRN);

    % Sample the sequence
    sampled_sequence = sample_sequence(caCode, Ts);

    % Calculate the shift index
    Tc = 1; % code period is 1 ms
    index = floor(rem_tot / Tc * length(sampled_sequence));

    % Shift the sequence
    % TODO not changing sequence 
    x_tau_seq = circshift(sampled_sequence, -index);
end

% function that accepts a time of travel, PRN, and sample rate it returns 
% the PRN sequence consistent with that state 
% function [earlyCode, lateCode, promptCode, remCodePhase, blksize] = get_x(settings, PRN, remCodePhase)
%     earlyLateSpc  = settings.dllCorrelatorSpacing;
%     codeFreq      = settings.codeFreqBasis;
%     codePhaseStep = codeFreq / settings.samplingFreq;
%     blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
%     % Get a vector with the C/A code sampled 1x/chip
%     caCode = generateCAcode(PRN);
%     % Then make it possible to do early and late versions
%     caCode = [caCode(1023) caCode caCode(1)];
% 
%     % Define index into early code vector
%     tcode       = (remCodePhase-earlyLateSpc) : ...
%                   codePhaseStep : ...
%                   ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
%     tcode2      = ceil(tcode) + 1;
%     earlyCode   = caCode(tcode2);
%     
%     % Define index into late code vector
%     tcode       = (remCodePhase+earlyLateSpc) : ...
%                   codePhaseStep : ...
%                   ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
%     tcode2      = ceil(tcode) + 1;
%     lateCode    = caCode(tcode2);
%     
%     % Define index into prompt code vector
%     tcode       = remCodePhase : ...
%                   codePhaseStep : ...
%                   ((blksize-1)*codePhaseStep+remCodePhase);
%     tcode2      = ceil(tcode) + 1;
%     promptCode  = caCode(tcode2);
%     
%     remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;
%   
% 
% end