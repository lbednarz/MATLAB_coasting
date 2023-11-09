classdef signalGenerator 
    properties
        lam
        f_IF
        Tacc
        c
    end

    methods 
        function obj = signalGenerator(settings)
            obj.lam = settings.lam;
            obj.f_IF = settings.f_IF;
            obj.Tacc = settings.Tacc;
            obj.c = settings.c;
        end 

        % based on: 
        %         gen_local_replica = @(x, range, I, T, dtr, dtsv, f_IF, isCosine) ...
        %     (x .* (isCosine .* cos(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc)) + ...
        %     (1-isCosine) .* sin(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc))));

        function localReplica = genLocalReplica(a, x, range, dtr, I, T, dtsv, isCosine)
            if isCosine == 1 
                localReplica = ...
                    a * x.*cos(2*pi/signalGenerator.lam * (range - I + T + signalGenerator.c*(dtr - dtsv) + ... 
                    signalGenerator.lam*signalGenerator.f_IF*signalGenerator.Tacc));
            else 
                localReplica = ... 
                    a * x.*sin(2*pi/signalGenerator.lam * (range - I + T + signalGenerator.c*(dtr - dtsv) + ...
                    signalGenerator.lam*signalGenerator.f_IF*signalGenerator.Tacc));
            end 
        end

        function RPrimeFunc = genRPrime(epsilon)
          RPrimeFunc = (epsilon < 0) .* 1/settings.codeFreqBasis + ...
            (epsilon > 0) .* -1/settings.codeFreqBasis;
        end

        function RFunc = genR(epsilon)
            RFunc = (abs(epsilon) < Tc) * ...
                         (1 - abs(epsilon)/Tc);
        end 

        function nominalSignal = genNominalSignal(a, epsilon, dtheta, isCosine)
            if isCosine == 1 
                nominalSignal = a * signalGenerator.grnR(epsilon) * cos(dtheta);
            else
                nominalSignal = a * signalGenerator.genR(epsilon) * sin(dtheta);
            end
        end

        function HFuncI = getHxI(a, epsilon, r_r, r_sv, dtheta)
            R = signalGenerator.genR(epsilon);
            Rp = signalGenerator.genRPrime(epsilon);
            dr = @(r_term, sv_term, r_r, r_sv)...
                 -2*(sv_term - r_term)/norm(r_sv-r_r) + (sv_term - r_term)/(norm(r_sv-r_r)^3);
            dterm = @(a, epsilon, dtheta, It) a/signalGenerator.c*Rp*cos(dtheta) ...
                - It*2*pi*a/signalGenerator.lam*R*sin(dtheta);

            HFuncI = ...
                [(a/signalGenerator.c * R_prime_func(epsilon) *cos(dtheta) ...
                    - 2*pi*a/signalGenerator.lam*R_func(epsilon)*sin(dtheta))*dr(r_r(1), r_sv(1), r_r, r_sv), ...
                 (a/signalGenerator.c * R_prime_func(epsilon) *cos(dtheta) ...
                    - 2*pi*a/signalGenerator.lam*R_func(epsilon)*sin(dtheta))*dr(r_r(2), r_sv(2), r_r, r_sv), ...
                 (a/signalGenerator.c * R_prime_func(epsilon) *cos(dtheta) ...
                    - 2*pi*a/signalGenerator.lam*R_func(epsilon)*sin(dtheta))*dr(r_r(3), r_sv(3), r_r, r_sv), ...
                 dterm(a, epsilon, dtheta, -1), ...
                 dterm(a, epsilon, dtheta,  1), ... 
                 dterm(a, epsilon, dtheta,  1), ... 
                 R_func(epsilon)*cos(dtheta)];
        end

        function HFuncQ = getHxQ(a, epsilon, r_r, r_sv, dtheta) 
            R = signalGenerator.genR(epsilon);
            Rp = signalGenerator.genRPrime(epsilon);
            dr = @(r_term, sv_term, r_r, r_sv)...
                 -2*(sv_term - r_term)/norm(r_sv-r_r) + (sv_term - r_term)/(norm(r_sv-r_r)^3);
            dterm = @(a, epsilon, dtheta, It) a/signalGenerator.c*Rp*sin(dtheta) ...
                 + It*2*pi*a/signalGenerator.lam*R*cos(dtheta);

            HFuncQ = ...
                [(a/signalGenerator.c * R_prime_func(epsilon) * sin(dtheta) ...
                    + 2*pi*a/signalGenerator.lam*R_func(epsilon) * cos(dtheta))*dr(r_r(1), r_sv(1), r_r, r_sv), ...
                 (a/signalGenerator.c * R_prime_func(epsilon) * cos(dtheta) ...
                    + 2*pi*a/signalGenerator.lam*R_func(epsilon) * sin(dtheta))*dr(r_r(2), r_sv(2), r_r, r_sv), ...
                 (a/signalGenerator.c * R_prime_func(epsilon) * cos(dtheta) ...
                    + 2*pi*a/signalGenerator.lam*R_func(epsilon) * sin(dtheta))*dr(r_r(3), r_sv(3), r_r, r_sv), ...
                 dterm(a, epsilon, dtheta, -1), ...
                 dterm(a, epsilon, dtheta,  1), ... 
                 dterm(a, epsilon, dtheta,  1), ... 
                 R_func(epsilon)*cos(dtheta)];
        end
        
        function CN0_dBhz = getCN0(Ip, Qp, Ie, Qe, Il, Ql) 
           CN0_dBhz= 10 * log10((Ip^2 + Qp^2)^2 / (0.5 * ((Ie^2 + Qe^2)^2 + ...
                    (Il^2 + Ql^2)^2)));
        end

    end

end


% gen_local_replica = @(x, range, I, T, dtr, dtsv, f_IF, isCosine) ...
%     (x .* (isCosine .* cos(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc)) + ...
%     (1-isCosine) .* sin(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc))));
% 
% % Define R(epsilon) 
% %samples_per_chip = .5*settings.samplingFreq*1/settings.codeFreqBasis; 
% % R_func = @(epsilon) (abs(epsilon) < Tc) * settings.codeLength * ...
% %                         .5*settings.samplingFreq/settings.codeFreqBasis * (1 - abs(epsilon)/Tc);
% R_func = @(epsilon) (abs(epsilon) < Tc) * ...
%                          (1 - abs(epsilon)/Tc);
% 
% 
% % Define R'(epsilon)
% % samples per chip is calculated as
% % sample_rate/(sample_rate/chip_frequency) = (sample/s) / (chips/s) =
% % chip_frequency
% R_prime_func = @(epsilon) (epsilon < 0) .* 1/settings.codeFreqBasis + ...
%                           (epsilon > 0) .* -1/settings.codeFreqBasis;
% 
% % define function for nominal signal measurement 
% gen_nominal_signal = @(a, epsilon, dtheta, isCosine) ...
%     a * R_func(epsilon) .* (isCosine .* cos(dtheta) + ...
%     (1-isCosine) .* sin(dtheta));
% 
% % shorthand for the derivative of of the function args wrt r_r
% dr = @(r_term, sv_term, r_r, r_sv)...
%      -1*(sv_term - r_term)/norm(r_sv-r_r);
% 
% % shorthand for the derivative of funtion wrt I, T, dtr 
% dterm = @(a, epsilon, dtheta, isCosine,It) ...
%          a/c*R_prime_func(epsilon)*(isCosine*cos(dtheta)+(1-isCosine)*(sin(dtheta))) ...
%             + It*2*pi*a/lam*R_func(epsilon)*(isCosine*sin(dtheta)+(1-isCosine)*(cos(dtheta)));
% 
% dtermP = @(a, dtheta, isCosine,It) ...
%       It*2*pi*a/lam*R_func(0)*(isCosine*sin(dtheta)+(1-isCosine)*(cos(dtheta)));
% 
% % here we are trying to enforce 1D, body frame movement along r_los 
% % this means dr_r2,r3 are zero, but a nominal range must be provided 
% % in the rotation matrix, r_r3 is along the los 
% % state vector order is [r_r, dtr, I, T, a];
% h_func_I = @(a, epsilon, r_r, r_sv, dtheta) ...
%     [a/c * R_prime_func(epsilon) * ...
%         dr(r_r(3), r_sv(3), r_r, r_sv)*cos(dtheta) ...
%             - 2*pi*a*R_func(epsilon)*sin(dtheta)*dr(r_r(3), r_sv(3), r_r, r_sv), ...
%      dterm(a, epsilon, dtheta, 1, -1), ...
%      dterm(a, epsilon, dtheta, 1,  1), ... 
%      dterm(a, epsilon, dtheta, 1,  1), ... 
%      R_func(epsilon)*cos(dtheta)];
% 
% h_func_Q = @(a, epsilon, r_r, r_sv, dtheta) ...
%     [a/c * R_prime_func(epsilon) * ...
%         dr(r_r(3), r_sv(3), r_r, r_sv)*sin(dtheta) ...
%             + 2*pi*a*R_func(epsilon)*cos(dtheta)*dr(r_r(3), r_sv(3), r_r, r_sv), ...
%      dterm(a, epsilon, dtheta, 0, -1), ...
%      dterm(a, epsilon, dtheta, 0,  1), ... 
%      dterm(a, epsilon, dtheta, 0,  1), ... 
%      R_func(epsilon)*sin(dtheta)];
% 
% h_func_P_I = @(a, r_r, r_sv, dtheta) ... % epsilon can only be zero
%     [2*pi*a/lam*R_func(0)*cos(dtheta)*dr(r_r(3), r_sv(3), r_r, r_sv), ...
%      dtermP(a, dtheta, 1, -1), ...
%      dtermP(a, dtheta, 1,  1), ... 
%      dtermP(a, dtheta, 1,  1), ... 
%      R_func(0)*sin(dtheta)];