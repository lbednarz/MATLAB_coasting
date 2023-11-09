classdef kfTracker
    properties 
        settings
        IMUct 
        signalGenerator
        xhat
        Phat 
        xbar 
        Pbar
        sv_pos
        dtsv
    end

    methods 
        function obj = kfTracker(settings_in, signalgenerator_in, priors)
            obj.settings = settings_in;
            obj.IMUct = getCtIMU(settings.tau_b_acc, settings.tau_b_gy);
            obj.signalGenerator = signalgenerator_in;
            obj.sv_pos = priors.sv_pos; 
            obj.dtsv = priors.dtsv;
        end

        function IMUct = getCtIMU(tau_b_acc, tau_b_gy)
            [Ac, Bc, Gc] = getIMUModelContinuous(tau_b_acc, tau_b_gy);
            % add iono, tropo, and signal amplitude must be appended 
            [an,am] = size(Ac);
            [~,bm] = size(Bc);
            [gn, gm] = size(Gc);
            Ac = [Ac zeros(an, 3);
                  zeros(1,am) -1/tau_I 0 0;
                  zeros(1,am) 0 -1/tau_T 0;
                  zeros(1,am) 0 0 1];
            Bc = [Bc; zeros(3,bm)];
            Gc = [Gc zeros(gn, 3); zeros(3, gm) eye(3)];
            syms phi theta psi_var phi_dot theta_dot psi_var_dot fb1 fb2 fb3 lat lon;
            IMUct = matlabFunction(Ac, Bc, Gc, 'Vars', ...
                                 [phi, theta, psi_var, ...
                                  phi_dot, theta_dot, psi_var_dot, ...
                                  fb1, fb2, fb3, ...
                                  lat, lon]);
        end

        function new_sys = time_update(x_lin, u)
            % get dynamics models 
            phi_k = x_lin(1);
            theta_k = x_lin(2);
            psi_k = x_lin(3);
            phi_dot_k = x_lin(4);
            theta_dot_k = x_lin(5);
            psi_dot_k = x_lin(6);
            fb1_k = x_lin(7);
            fb2_k = x_lin(8);
            fb3_k = x_lin(9);

            lat = kfTracker.settings.lat;
            lon = kfTracker.settings.lon;

            % 15 state IMU 
            [Ac, Bc, Gc] = kfTracker.IMUct(phi_k, theta_k, psi_k, ...
                                  phi_dot_k, theta_dot_k, psi_dot_k, ...
                                  fb1_k, fb2_k, fb3_k, ... 
                                  lat, lon);
            C = [];
            D = 0;
            [~,am] = size(Ac);
            [~,bm] = size(Bc);
            [gn,gm] = size(Gc);
            Bc = [Bc; zeros(5,bm)]; % 5 comes from 2 state clock model, Iono, Tropo, and signal amplitude equations 
            Gc = [Gc zeros(gn,5);
                  zeros(5,gm) eye(5)];
            sysc = ss(Ac, Bc, C, D);
            sysd = c2d(sysc, kfTracker.settings.TsIMU, 'zoh');
            A = sysd.A;
            % non-IMU models (2 state clock, Iono, tropo, amplitude)
            A_other = [zeros(1,am) 1-kfTracker.settings.TsIMU/kfTracker.settings.tau_I zeros(1,4);
                       zeros(1,am+1) 1-kfTracker.settings.TsIMU/kfTracker.settings.tau_T zeros(1,3);
                       zeros(1,am+2) 1 kfTracker.settings.TsIMU 0;
                       zeros(1,am+3) kfTracker.settings.TsIMU 0;
                       zeros(1,am+4) 1];
            A = [A;A_other];
            B = sysd.B;  
            Q = Gc*kfTracker.settings.W*Gc.' * kfTracker.settings.TsIMU;
            
            kfTracker.xbar = A*kfTracker.xhat + B*u;
            kfTracker.Pbar = A*kfTracker.Phat*A.' + Q;  
            
            new_sys.xbar = kfTracker.xbar;
            new_sys.Pbar = kfTracker.Pbar;

        end

        function new_sys = msmtUpdate(rawSignal, rem_tau, PRN, lp_filter) 
            
            range = norm(kfTracker.sv_pos - kfTracker.xbar);
            dsv = kfTracker.dtsv;
            xk = kfTracker.xbar;

            % Jacobian and state scaling are performed 
            Dx = ones(size(A));
            Dxk(1) = 1/x_hat_LOS(3);
            Dy = 1/num_samples*ones(4,size(xk,1)); % num measurements is hardcoded here! 
            epsilon = kfTracker.settings.epsilon;

            early_code = get_xk(rem_tau-epsilon, PRN, ...
                    kfTracker.settings.sample_rate);
            late_code  = get_xk(rem_tau+epsilon*1000, PRN, ...
                    kfTracker.settings.sample_rate);
            prompt_code = get_xk(rem_tau, PRN,  ... 
                    kfTracker.settings.sample_rate);
            
            sg = kfTracker.signalGenerator; 

            LR_early_I  = sg.genLocalReplica(xk(20), early_code, range, ...
                xk(16), xk(18), xk(19), dsv, 1);
            LR_early_Q  = sg.genLocalReplica(xk(20), early_code, range, ... 
                xk(16), xk(18), xk(19), dsv, 0);
            LR_late_I   = sg.genLocalReplica(xk(20), late_code, range, ...
                xk(16), xk(18), xk(19), dsv, 1);
            LR_late_Q   = sg.genLocalReplica(xk(20), late_code, range, ... 
                xk(16), xk(18), xk(19), dsv, 0);
            LR_prompt_I = sg.genLocalReplica(xk(20), prompt_code, range, ...
                xk(16), xk(18), xk(19), dsv, 1);
            LR_prompt_Q = sg.genLocalReplica(xk(20), prompt_code, range, ...
                xk(16), xk(18), xk(19), dsv, 0); 
        
            % mix signals  
            I_E = LR_early_I  .* real(rawSignal);
            Q_E = LR_early_Q  .* imag(rawSignal);
            I_L = LR_late_I   .* real(rawSignal);
            Q_L = LR_late_Q   .* imag(rawSignal);
            I_P = LR_prompt_I .* real(rawSignal);
            Q_P = LR_prompt_Q .* imag(rawSignal);
            
            % low pass 
            I_E_bp = filter(lp_filter, I_E);
            Q_E_bp = filter(lp_filter, Q_E);
            I_L_bp = filter(lp_filter, I_L);
            Q_L_bp = filter(lp_filter, Q_L);
            I_P_bp = filter(lp_filter, I_P);
            Q_P_bp = filter(lp_filter, Q_P);
        
            % correlate
            I_E_filt = sum(I_E_bp);
            Q_E_filt = sum(Q_E_bp);
            I_L_filt = sum(I_L_bp);
            Q_L_filt = sum(Q_L_bp);
            I_P_filt = sum(I_P_bp);
            Q_P_filt = sum(Q_P_bp);
            
            % concatenate into measurement
            y_filt   = 1/num_samples * [I_E_filt;Q_E_filt;I_L_filt;Q_L_filt];
            % estimate measurement noise
            CN0dBHz = CN0_dBHz(I_P_filt, Q_P_filt, I_E_filt, Q_E_filt, I_L_filt, Q_L_filt);
            CN0 = 10^(CN0dBHz/10); % Carrier to noise density ratio
            R = diag(1/CN0*ones(4,1)); % Variance of measurement noise on 4 measurements 
            R = (inv(Dy))'*R*(inv(Dy)); % scaling the jacobian requires this

            HxIE = sg.getHxI(xk(20), -epsilon, xk(1:3), kfTracker.sv_pos, ...
                atan2(Q_E_filt,I_E_filt));
            HxIL = sg.getHxI(xk(20), epsilon, xk(1:3), kfTracker.sv_pos, ...
                atan2(Q_L_filt,I_L_filt));
            HxQE = sg.getHxQ(xk(20), -epsilon, xk(1:3), kfTracker.sv_pos, ...
                atan2(Q_E_filt,I_E_filt));
            HxQL = sg.getHxQ(xk(20), epsilon, xk(1:3), kfTracker.sv_pos, ...
                atan2(Q_L_filt,I_L_filt));

            % state vector for HxIE, for instance, is 
            % rx, ry, rz, I, T, dtr, a
            HstarTmp = [HxIE; HxQE; HxIL; HxQL];
            x_tmp = [xk(1:3), xk(16), xk(18:end)].';
            y_star = HstarTmp*x_tmp; 
            h_star = zeros(length(y_filt), length(xk));
            h_star(:,1:3) = HstarTmp(:,1:3);
            h_star(:,16:18) = HstarTmp(:,4:6);
            h_star(:,end) = HstarTmp(:,end);
            
            % generate nominal signal to generate linearized observtion
            ynom_I_E = sg.genNominalSignal(xk(end), -1*epsilon, atan2(Q_E_filt,I_E_filt), 1);
            ynom_Q_E = sg.genNominalSignal(xk(end), -1*epsilon, atan2(Q_E_filt,I_E_filt), 0);
            ynom_I_L = gen_nominal_signal(xk(end), epsilon, atan2(Q_L_filt,I_L_filt), 1);
            ynom_Q_L = gen_nominal_signal(xk(end), epsilon, atan2(Q_L_filt,I_L_filt), 0);
            
            % concatenate into nominal measurement 
            y_nom    = [ynom_I_E;ynom_Q_E;ynom_I_L;ynom_Q_L];
        
            % form adjusted measurement 
            y_star = Dy^-1*(y_filt - y_nom) + y_star; % jacobian scaling changes measurement 
        
%             % add zeros to accomodate this
%     h_star = [h_star_tmp(:,1), zeros(4,1), h_star_tmp(:,2), zeros(4,1), ...
%                 h_star_tmp(:,3:4), zeros(4,1), h_star_tmp(:,5)];    
    
            % run measurement update
            K = Pk * h_star' * inv(h_star * Pk * h_star' + R);
            innov = y_star - h_star*xk;
            xhatk = xk + K*(innov);
            P_ch = (eye(size(Pk)) - K*h_star);
            Phatk = P_ch*Pk;
        end

    end
end