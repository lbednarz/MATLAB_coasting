clear
close all

c = 2.99792458e8;
L5 = 1176.45e6;
lam5 = c/L5;

CN0dBHz = -0;
CN0 = 10.^(CN0dBHz/10);
R = 1/CN0;

%% CLOCK -- LPFRS

h0 = 1.5e-22;
h0 = h0*20;     % *20 w/ vibration
h_2 = 8.5e-32;

Sf = h0/2*L5^2;
Sg = 2*pi^2*h_2*L5^2;

Fc = [0 1; 0 0];
Gc = eye(2);
Qc = diag([Sf Sg]);
Hc = [1 0];

sysc = ss(Fc,Gc,Hc,[]);
[KESTc,Lc,Pc] = kalman(sysc,Qc,R,[]);
KFc = [1 0 0]*KESTc;
tf(KFc)
figure(1);
bode(KFc); grid

% xlabel('time [ms]', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('$\hat{P_{r_r}}$', 'Interpreter', 'latex', 'FontSize', 16);
% title('Plot of Receiver Postion Covariance vs. time', 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 14);

% Save as a PNG with 300 dots per inch (DPI)
print('bode', '-dpng', '-r300');

sig_phase_clk_deg = sqrt(Pc(1,1))*180/pi
sig_phase_clk_cm = sig_phase_clk_deg*25.5/360

CN0dBHz_range = -20:30;
plot_CN0_vs_Pc(sysc, Qc, CN0dBHz_range)




%% ACCELEROMETER -- Ellipse2

% g = 9.81; % m/s^2
% 
% VRW = 0.033; % m/s/sqrt(hr)
% Qaw = 100*(VRW/60)^2;  % (m/s^2)^2/Hz
% 
% tauab = 500; % sec
% sigab = 0.014; % mg
% sab = sigab/1000*g; % m/^2 
% Qab = 2*tauab*sab^2;
% 
% % states are phase, dr, dv, b
% % PHASE = f_o*[0 0, 1/c, 0] + wa - vsv/c
% Fi = [0 0 1/c 0; 0 1 0 0; 0 1 0 1;0 0 0 -1/tauab];
% Gi = [1 0 0;0 0 0; 0 1 0;0 0 1/tauab];
% Qi = diag([Sf Qaw Qab]);
% %Qi = Gi * Qi * Gi';
% %Hi = [1 0 0 0];
% 
% sysi = ss(Fi,Gi,Hi,[]);
% [KESTi,Li,Pi] = kalman(sysi,Qi,R,[]);
% KFi = [1 0 0 0]*KESTi;
% tf(KFi)
% figure(2); 
% %bode(KFi); grid
% 
% sig_phase_imu_deg = sqrt(Pi(1,1))*180/pi
% sig_phase_imu_cm = sig_phase_imu_deg*lam5/360
% 
% % construct range 
% %dr(k) = dr(k-1) + dv(k-1);
% CN0dBHz_range = -10:35;
% plot_CN0_vs_Pc(sysi, Qi, CN0dBHz_range)

% %% build priors 
% % close all; clear; clc; 
% addpath 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GNSS_SDR-master'
% addpath data
% 
% load("data\extras97.mat"); % should set a variable called "extra_out"
% load("data\trackingResults96.mat")
% channel_to_track = 8;                % channel we want to track
% 
% % initials for position and signal quality
% settings = initSettings();
% c = 2.99792458e8;             % speed of light in m/s
% L1 = 1575.42e6;               % L1 carrier frequency
% lam = c/L1;                   % carrier wavelength
% Tc = 1/1.023e6;               % chip duration
% Ts = settings.samplingFreq/2; % Hz sample rate
% Tacc = 1e-3;                  % accumulation period (seconds) 
% f_IF = settings.IF;           % intermediate frequency (Hz)
% 
% % intial estimates 
% % [x_hat, y_hat, z_hat, cdtr_hat] = priors_1.state; 
% priors_1 = extra_out.extras_epoch.epoch1;
% x_hat = priors_1.state;
% 
% % nominal values for local replica
% range = norm(priors_1.satpos(:,channel_to_track)-x_hat(1:3));
% tau = priors_1.tot*1000; % travel time in ms 
% rem_tau = tau - floor(tau); % retrieve just the time into current period
% epsilon = .5*Tc;
% dtsv = priors_1.dtsv(channel_to_track);
% 
% % other params for kalman filter 
% % signal power 
% PdBHz = -130;      % Signal power dB/Hz
% P = 10^(PdBHz/10); % convert to V? 
% %sig_P = 5e-7;      % V?
% sig_P = 50;      % V?  
% 
% 
% % Iono and tropo
% tau_I = 30*60;    % About a 30-minute time constant
% tau_T = 2*3600;   % About a 2-hour time constant
% sig_I = 2;        % Should be about between 1-5 meters
% sig_T = .5*sig_I; % Should be roughly 1/2 of Iono
% 
% % Accelerometer with bias
% sig_a = 1e5*57e-6; % SD accelerometer noise m/s^1.5
% sig_b = 1e3*14e-6; % Continuous SD of accelerometer bias RW noise
% tau_b = 3600;  % Time constant of accelerometer bias RW
% 
% % Clock TXCO
% sig_d = 1e4*5e-8 * sqrt(Ts); % Discrete clock drift noise variance
% 
% % Assuming values for Ts, tau_I, tau_T, tau_b, and g_adj are provided:
% dt = 1; % time between integration periods 
% A = [1, dt, 0,  0,       0,        0,            0,    0;
%      0,  1, 0,  0,       0,        0,          -dt,    0;
%      0,  0, 1, dt,       0,        0,            0,    0;
%      0,  0, 0,  1,       0,        0,            0,    0;
%      0,  0, 0,  0,   1-dt/tau_I,   0,            0,    0;
%      0,  0, 0,  0,       0,     1-dt/tau_T,      0,    0;
%      0,  0, 0,  0,       0,        0,      1-dt/tau_b, 0;
%      0,  0, 0,  0,       0,        0,            0,    1];
% 
% % control inputs are r_dot_nominal, v_nominal, v_dot_nominal, and gadj 
% B = [dt, -dt,  0,   0;
%       0,   0, dt, -dt;
%       0,   0,  0,   0;
%       0,   0,  0,   0;
%       0,   0,  0,   0;
%       0,   0,  0,   0;
%       0,   0,  0,   0;
%       0,   0,  0,   0];
% 
% G = diag([-dt, 0.5*dt^2, 1, dt, dt, dt, dt]);
% [~,m] = size(G);
% G = vertcat(zeros(1,m),G);
% W = diag([sig_a^2, sig_d^2, sig_d^2, sig_I^2, sig_T^2, sig_b^2, sig_P^2]);
% Q = G * W * G';
% Q(1,1)=1;
% 
% H = eye(size(A));
% 
% sysc = ss(A,B,H,[]);
% 
% size(sysc.A)
% size(Q)
% 
% CN0dBHz_range = -10:35;
% plot_CN0_vs_Pc(sysc, Q, CN0dBHz_range)