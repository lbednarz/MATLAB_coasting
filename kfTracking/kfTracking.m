%% load in prior infomation from PLL/DLL
% for this demo, we are interested in channels 2,3,8,9
close all; clear; clc; 
addpath 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GNSS_SDR-master'
addpath data

load("data\extras97.mat"); % should set a variable called "extra_out"
load("data\trackingResults96.mat")
dataAdaptCoeff = 1; 
channel_to_track = 8;                % channel we want to track
PRN = 13;                            % PRN in channel
frame_starts = extra_out.sample_num; % flag for each subframe 
byte_frame_start = trackResults(channel_to_track).absoluteSample(frame_starts(1,channel_to_track)); % should be 220614790
settings = initSettings();

% get RF file to track 
fileName = 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GPSStatic_BW8_gain65_10Msps__REAL.bin';
[fid, message] = fopen(fileName, 'rb');

% Move to location of first prior 
  fseek(fid, ...
    dataAdaptCoeff*byte_frame_start, ...
    'bof');

% parameters for reading raw data 
remCodePhase = 0;
codeFreq      = settings.codeFreqBasis;
codePhaseStep = codeFreq / settings.samplingFreq;

%% build priors 
% initials for position and signal quality
c = 2.99792458e8;             % speed of light in m/s
L1 = 1575.42e6;               % L1 carrier frequency
lam = c/L1;                   % carrier wavelength
Tc = 1/1.023e6;               % chip duration
Ts = settings.samplingFreq/2; % Hz sample rate
Tacc = 1e-3;                  % accumulation period (seconds) 
f_IF = settings.IF;           % intermediate frequency (Hz)

% intial estimates 
% [x_hat, y_hat, z_hat, cdtr_hat] = priors_1.state; 
priors_1 = extra_out.extras_epoch.epoch1;
x_hat = priors_1.state;

% nominal values for local replica
range = norm(priors_1.satpos(:,channel_to_track)-x_hat(1:3));
tau = priors_1.tot*1000; % travel time in ms 
rem_tau = tau - floor(tau); % retrieve just the time into current period
epsilon = .5*Tc;
dtsv = priors_1.dtsv(channel_to_track);

% other params for kalman filter 
% signal power 
PdBHz = -130;      % Signal power dB/Hz
P = 10^(PdBHz/10); % convert to V? 
sig_P = 5e-7;      % V?  

% Iono and tropo
tau_I = 30*60;    % About a 30-minute time constant
tau_T = 2*3600;   % About a 2-hour time constant
sig_I = 2;        % Should be about between 1-5 meters
sig_T = .5*sig_I; % Should be roughly 1/2 of Iono

% Accelerometer with bias
sig_a = 57e-6; % SD accelerometer noise m/s^1.5
sig_b = 14e-6; % Continuous SD of accelerometer bias RW noise
tau_b = 3600;  % Time constant of accelerometer bias RW

% Clock TXCO
sig_d = 5e-8 * sqrt(Ts); % Discrete clock drift noise variance

% convert ECEF pos to LOS basis 
R_EB = LOSRotationMatrix(x_hat(1:3), priors_1.satpos(:,channel_to_track));
x_hat_LOS = R_EB*x_hat(1:3); % x_hat_LOS(1 and 2) are constant 

% initial covairances
% take state vector r_r, v_r, dtr, dtrdot, I, T, b, P
P0 = priors_1.Q; % assumes states are [x, y, z, c*dtr]
P0(1:3,1:3) = R_EB*P0(1:3,1:3); % rotate into body frame
Prdt = [P0(2,2), P0(2,4); P0(4,2), P0(4,4)].*1000; % 1D pos and clock uncertainty
Pv0 = 1; Pdtrdot0 = 5e-5; PT0 = 3^2; PI0 = 3^2; Pb0 = 14e-3^2; PP0 = 10^2; % TODO PP0 might be crazy 
P0 = diag([Prdt(1,1), Pv0, Prdt(2,2), Pdtrdot0, PI0, PT0, Pb0, PP0]);
P0(1,3) = Prdt(1,2); P0(3,1) = Prdt(1,2);

% initial states 
rk = x_hat_LOS(3);
vk = 0;
dtrk = x_hat(4)/c;   % Initial receiver clock bias
dtrdotk = 2*sig_d;   % Initial receiver clock drift
Ik = 0;
Tk = priors_1.tropo;
bk = 14e-6;          % Initial guess at bias
ak = sqrt(2*P);

%% measurement models 
% get spreading code for specified sample rate, PRN, and time of travel
% using get_x(rem_tot, PRN, Ts)
% x comes from get_x function
% isCosine allows for making I and Q replicas from 1 function 
gen_local_replica = @(x, range, I, T, dtr, dtsv, f_IF, isCosine) ...
    x .* (isCosine .* cos(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc)) + ...
    (1-isCosine) .* sin(2*pi/lam * (range - I + T + c*(dtr - dtsv) + lam*f_IF*Tacc)));

% Define R(epsilon) 
R_func = @(epsilon) (abs(epsilon) < Tc) .* (1 - abs(epsilon)/Tc);

% Define R'(epsilon)
R_prime_func = @(epsilon) (epsilon < 0) .* 1.023e6 + ...
                          (epsilon > 0) .* -1.023e6;

% define function for nominal signal measurement 
gen_nominal_signal = @(a, epsilon, range, I, T, c, dtr, dtsv, isCosine) ...
    a * R_func(epsilon) .* (isCosine .* cos(2*pi/lambda * (range - I + T + c*(dtr - dtsv))) + ...
    (1-isCosine) .* sin(2*pi/lambda * (range - I + T + c*(dtr - dtsv))));

% here we are trying to enforce 1D, body frame movement along r_los 
% this means dr_r2,r3 are zero, but a nominal range must be provided 
% in the rotation matrix, r_r3 is along the los 
h_func_I = @(a, epsilon, r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3) ...
    [R_func(epsilon), ... 
    -a/c * R_prime_func(epsilon) * (r_sv3 - r_r3)/norm([r_sv1;r_sv2;r_sv3]-[r_r1;r_r2;r_r3]), ...
    a * R_prime_func(epsilon) / c, ... 
    a * R_prime_func(epsilon) / c, ... 
    a * R_prime_func(epsilon) / c];

h_func_Q = @(a, epsilon, r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3) ...
    [0, ... 
    -2*pi*a/lam * R_func(epsilon) * (r_sv3 - r_r3)/norm([r_sv1;r_sv2;r_sv3]-[r_r1;r_r2;r_r3]), ...
    -2*pi*a/lam * R_func(epsilon), ... 
    2*pi*a/lam * R_func(epsilon), ... 
    2*pi*a/lam * R_func(epsilon)];

% CN0 estimator - only thing using I_p for now
CN0_dBHz = @(Ip, Qp, Ie, Qe, Il, Ql) ...
    10 * log10((Ip^2 + Qp^2) / (0.5 * ((Ie^2 + Qe^2) + ...
    (Il^2 + Ql^2)) - (Ip^2 + Qp^2))); 

% low pass filter for after signal accumulation 
fs = 10e6; % Sampling frequency, for example 10 MHz
fpass = 2.5e6; % Passband frequency
fstop = 3e6; % Stopband frequency
atten_stopband = 60; % Attenuation in the stopband in dB
% Apass = 1;     % Passband ripple (in dB) and add 'PassbandRipple', Apass,
% if you want this later 

% Design a low-pass FIR filter
lp_filter = designfilt('lowpassfir', 'PassbandFrequency', fpass, ...
                       'StopbandFrequency', fstop, 'StopbandAttenuation', atten_stopband, ...
                       'SampleRate', fs, 'DesignMethod', 'kaiserwin');


%% dynamic models 
% take state vector r_r, v_r, dtr, dtrdot, I, T, b, P
% r_{k} = r_{k-1} + (v_{k-1} + (r_dot_nominal_{k-1} - v_nominal_{k-1}))*Ts
% v_{k} = v_{k-1} + (-1*b + (v_dot_nominal_{k-1}-g_adj) - normrnd(0, sig_a))*Ts
% dtr_{k} = dtr_{k-1} + Ts*dtr_dot_{k-1} + 0.5*normrnd(0, sig_d)*Ts^2
% dtr_dot_{k} = dtr_dot_{k-1} + normrnd(0, sig_d)
% I_{k} = I_{k-1} + (-1/tau_I * I_{k} + normrnd(0, sig_I))*Ts
% T_{k} = T_{k-1} + (-1/tau_T * T_{k} + normrnd(0, sig_T))*Ts
% b_{k} = b_{k-1} + (-1/tau_b * b_{k-1} + normrnd(0, sig_b))*Ts
% P_{k} = P_{k-1} + normrnd(0, sig_P)*Ts

% Assuming values for Ts, tau_I, tau_T, tau_b, and g_adj are provided:

A = [1, Ts, 0,  0,       0,        0,            0,    0;
     0,  1, 0,  0,       0,        0,          -Ts,    0;
     0,  0, 1, Ts,       0,        0,            0,    0;
     0,  0, 0,  1,       0,        0,            0,    0;
     0,  0, 0,  0,   1-Ts/tau_I,   0,            0,    0;
     0,  0, 0,  0,       0,     1-Ts/tau_T,      0,    0;
     0,  0, 0,  0,       0,        0,      1-Ts/tau_b, 0;
     0,  0, 0,  0,       0,        0,            0,    1];

% control inputs are r_dot_nominal, v_nominal, v_dot_nominal, and gadj 
B = [Ts, -Ts,  0,   0;
      0,   0, Ts, -Ts;
      0,   0,  0,   0;
      0,   0,  0,   0;
      0,   0,  0,   0;
      0,   0,  0,   0;
      0,   0,  0,   0;
      0,   0,  0,   0];

G = diag([-Ts, 0.5*Ts^2, 1, Ts, Ts, Ts, Ts]);
[~,m] = size(G);
G = vertcat(zeros(1,m),G);
W = diag([sig_a^2, sig_d^2, sig_d^2, sig_I^2, sig_T^2, sig_b^2, sig_P^2]);
Q = G * W * G';

% need to implement a CN0 estimator here  
% currently, I will have to implement a prompt correlator to estimate the
% signal power. The early and late correlators will be used to estimate the
% noise power. The final scheme will be: 
% P_s = norm([I_p; Q_p])^2
% P_e = norm([I_e; Q_e])^2
% P_l = norm([I_l;Q_l])^2
% P_n = (P_e + P_l)/2
% CN0dBHz = 10*log10(Ps/Pn)
% CN0 = 10^(CN0dBHz/10); % Carrier to noise density ratio
% R = 1/CN0; % Variance of measurement noise
% sig_y = sqrt(R); % SD of measurement noise
% Anonymous function to compute CN0 in dBHz
% r_los = @(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3) ...
%     ([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]) / ...
%     norm([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]);

%% Tracking loop

% get navbits 
nb = extra_out.extra_nb{1,channel_to_track};
nb(nb == 0) = -1; % encode NB into 1s and -1s
sim_time = 30; % simulation time (s)

% get first block of raw data
blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
[rawSignal, samplesRead] = fread(fid, ...
    dataAdaptCoeff*blksize, settings.dataType);

rawSignal = rawSignal'; 

% raw signal w/o navbit
% repeat for each measurement epoch. Next batch will get a new navbit that
% will be applicable for 20 code periods.
rawSignal = nb(1)*rawSignal;

% % nominal values for local replica
% range = norm(priors_1.satpos(:,channel_to_track)-x_hat(1:3));
% tau = priors_1.tot*1000; % travel time in ms 
% rem_tau = tau - floor(tau); % retrieve just the time into current period
% epsilon = .5*Tc;
% dtsv = priors_1.dtsv(channel_to_track);
% % initial states 
% r0 = x_hat_LOS(3);
% v0 = 0;
% dtr0 = x_hat(4)/c;   % Initial receiver clock bias
% dtrdot0 = 2*sig_d;   % Initial receiver clock drift
% I0 = 0;
% T0 = priors_1.tropo;
% b0 = 14e-6;          % Initial guess at bias
% a0 = sqrt(2*P);
Pk = P0;

% for i = 1:sim_time*1000 - 1 % each i is 1 ms  
    % TODO: may have to incorporate remainder of code/carrier phase
    % make correlators 
    early_code  = get_x(rem_tot-epsilon, PRN, Ts);
    late_code   = get_x(rem_tot+epsilon, PRN, Ts);
    prompt_code = get_x(rem_tot, PRN, Ts);
    LR_early_I  = gen_local_replica(early_code, range, Ik, Tk, dtrk, dtsv, f_IF, 1);
    LR_early_Q  = gen_local_replica(early_code, range, Ik, Tk, dtrk, dtsv, f_IF, 0);
    LR_late_I   = gen_local_replica(late_code, range, Ik, Tk, dtrk, dtsv, f_IF, 1);
    LR_late_Q   = gen_local_replica(late_code, range, Ik, Tk, dtrk, dtsv, f_IF, 0);
    LR_prompt_I = gen_local_replica(prompt_code, range, Ik, Tk, dtrk, dtsv, f_IF, 1);
    LR_prompt_Q = gen_local_replica(prompt_code, range, Ik, Tk, dtrk, dtsv, f_IF, 0);

    % correlate 
    I_E = sum(LR_early_I  .* rawSignal);
    Q_E = sum(LR_early_Q  .* rawSignal);
    I_L = sum(LR_late_I   .* rawSignal);
    Q_L = sum(LR_late_Q   .* rawSignal);
    I_P = sum(LR_prompt_I .* rawSignal);
    Q_P = sum(LR_prompt_Q .* rawSignal);
    
    % low pass 
    I_E_filt = filtfilt(lp_filter, I_E);
    Q_E_filt = filtfilt(lp_filter, Q_E);
    I_L_filt = filtfilt(lp_filter, I_L);
    Q_L_filt = filtfilt(lp_filter, Q_L);
    I_P_filt = filtfilt(lp_filter, I_P);
    Q_P_filt = filtfilt(lp_filter, Q_P);

    % estimate measurement noise
    CN0dBHz = CN0_dBHz(I_P_filt, Q_P_filt, I_E_filt, Q_E_filt, I_L_filt, Q_L_filt);
    CN0 = 10^(CN0dBHz/10); % Carrier to noise density ratio
    R = 1/CN0; % Variance of measurement noise

    % build linearized observation matrix for Kalman filter
    


