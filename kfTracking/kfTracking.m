%% load in prior infomation from PLL/DLL
% for this demo, we are interested in channels 2,3,8,9
close all; clear; clc; beep off;
%addpath 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GNSS_SDR-master'
addpath data

load("data\extras97.mat"); % should set a variable called "extra_out"
load("data\trackingResults96.mat")
dataAdaptCoeff = 1; 
channel_to_track = 8;                % channel we want to track
PRN = 13;                            % PRN in channel
frame_starts = extra_out.sample_num; % flag for each subframe 
byte_frame_start = trackResults(channel_to_track).absoluteSample(frame_starts(1,channel_to_track)); % should be 220614790
settings = set_settings();

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
c = settings.c;               % speed of light in m/s
L1 = settings.carrier_freq;   % signal carrier frequency
lam = c/L1;                   % carrier wavelength
Tc = settings.codeFreqBasis;  % chip duration
Ts = settings.samplingFreq/2; % Hz sample rate
Tacc = settings.Tacc;         % accumulation period (seconds) 
f_IF = settings.IF;           % intermediate frequency (Hz)
num_samples = .5*settings.samplingFreq*1e-3;

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
sig_P = 50;      % V?  

% convert ECEF pos to LOS basis 
R_EB = LOSRotationMatrix(x_hat(1:3), priors_1.satpos(:,channel_to_track))';
x_hat_LOS = R_EB*x_hat(1:3); % x_hat_LOS(1 and 2) are constant 
sv_pos_LOS = R_EB*priors_1.satpos(:,channel_to_track); % SV pos in new coordinate frame

% initial covairances
% take state vector r_r, v_r, dtr, dtrdot, I, T, b, P
P0 = priors_1.Q; % assumes states are [x, y, z, c*dtr]
P0(1:3,1:3) = R_EB*P0(1:3,1:3); % rotate into body frame
Prdt = [P0(2,2)*1000, P0(2,4); P0(4,2), P0(4,4)]; % 1D pos and clock uncertainty
Pv0 = .5^2; Pdtrdot0 = 5e-5; PT0 = 3^2; PI0 = 3^2; Pb0 = (14e-3)^2; %PP0 = 1e-7^2;
PP0 = 100^2; %signal initial amplitude - anybody's guess from ADC  
P0 = diag([Prdt(1,1), Pv0, Prdt(2,2)/1000, Pdtrdot0, PI0, PT0, Pb0, PP0]);
P0(1,3) = Prdt(1,2); P0(3,1) = Prdt(1,2);

% initial states 
rk = x_hat_LOS(3);
vk = 0;
dtrk = x_hat(4)/c;   % Initial receiver clock bias
dtrdotk = 2*sig_d;   % Initial receiver clock drift
Ik = 0;
Tk = priors_1.tropo;
bk = 14e-6;          % Initial guess at bias
%ak = sqrt(2*P);
ak = 10; 

% rotate gravitational impact into body frame 
% Given Latitude in DMS format
lat_deg = 40;
lat_min = 48;
lat_sec = 33.2064;

% Convert Latitude from DMS to Decimal Degrees
lat_decimal = lat_deg + lat_min/60 + lat_sec/3600;

% Given Longitude in DMS format
lon_deg = 113; % Taking magnitude first, we'll apply the sign later
lon_min = 36;
lon_sec = 11.4036;

% Convert Longitude from DMS to Decimal Degrees
lon_decimal = -(lon_deg + lon_min/60 + lon_sec/3600); % Negative because it's west

lat_radians = deg2rad(lat_decimal);
lon_radians = deg2rad(lon_decimal);

R_ENU_to_ECEF = @(phi, lambda) ... % phi=lat, lambda = lon
    [-sin(lambda),           -sin(phi)*cos(lambda),    cos(phi)*cos(lambda);
      cos(lambda),           -sin(phi)*sin(lambda),    cos(phi)*sin(lambda);
      0,                      cos(phi),                 sin(phi)];

R_UE = R_ENU_to_ECEF(lat_radians,lon_radians);
g_ENU = [0;0;-9.81];
g_ECEF = R_UE*g_ENU;
g_adj = R_EB*g_ECEF; 
g_adj = g_adj(3); % grabbing part along LOS

%% measurement models 

% CN0 estimator - only thing using I_p for now
CN0_dBHz = @(Ip, Qp, Ie, Qe, Il, Ql) ...
    10 * log10((Ip^2 + Qp^2)^2 / (0.5 * ((Ie ^2 + Qe^2)^2 + ...
    (Il^2 + Ql^2)^2))); 

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

%% Tracking loop

% TODO reconsider how to grab data 


% grab I and Q components 
rawSignal = rawSignal';
rawSignal1=rawSignal(1:2:end);
rawSignal2=rawSignal(2:2:end);
rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector

% set initial covariance and states
Pk = P0;
xk = [rk, vk, dtrk, dtrdotk, Ik, Tk, bk, ak]';
counter = 1;
Pk_r = [sqrt(P0(1,1))];

% Jacobian and state scaling are performed 
Dx = ones(size(A));
Dx(1) = 1/x_hat_LOS(3);
Dy = 1/num_samples*ones(4,size(xk,1)); % num measurements is hardcoded here! 

% adjust dynamic models (measurement models are adjusted in the loop)  
B = Dx^-1*B;
A = Dx^-1*A*Dx;
Pk = Dx^-1*Pk*(Dx^-1)';
Q = Dx^-1*Q*(Dx^-1)';


N = 1000; % Number of points
data = zeros(1, N); % Storage for sequence
cn0_list = [];
range_l = [];


for i = 1:1000 % each i is 1 ms  
    % TODO: may have to incorporate remainder of code/carrier phase
    % make correlators 
    early_code  = get_x(rem_tau-epsilon*1000, PRN, Ts);
    late_code   = get_x(rem_tau+epsilon*1000, PRN, Ts);
    prompt_code = get_x(rem_tau, PRN, Ts);
    LR_early_I  = gen_local_replica(early_code, range, xk(5), xk(6), xk(3), ...
                                        dtsv, f_IF, 1);
    LR_early_Q  = gen_local_replica(early_code, range, xk(5), xk(6), xk(3), ... 
                                        dtsv, f_IF, 0);
    LR_late_I   = gen_local_replica(late_code, range, xk(5), xk(6), xk(3), ...
                                        dtsv, f_IF, 1);
    LR_late_Q   = gen_local_replica(late_code, range, xk(5), xk(6), xk(3), ...
                                        dtsv, f_IF, 0);
    LR_prompt_I = gen_local_replica(prompt_code, range, xk(5), xk(6), xk(3), ...
                                        dtsv, f_IF, 1);
    LR_prompt_Q = gen_local_replica(prompt_code, range, xk(5), xk(6), xk(3), ...
                                        dtsv, f_IF, 0);

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
    % TODO - is the still fair since I'm scaling the jacobian? 
    y_filt   = 1/num_samples * [I_E_filt;Q_E_filt;I_L_filt;Q_L_filt];

    % estimate measurement noise
    CN0dBHz = CN0_dBHz(I_P_filt, Q_P_filt, I_E_filt, Q_E_filt, I_L_filt, Q_L_filt);
    CN0 = 10^(CN0dBHz/10); % Carrier to noise density ratio
    R = diag(1/CN0*ones(4,1)); % Variance of measurement noise on 4 measurements 
    R = (inv(Dy))'*R*(inv(Dy)); % scaling the jacobian requires this 
    cn0_list = [cn0_list, CN0dBHz];

    % build linearized observation matrix for Kalman filter
    a = xk(8);
    h_I_E = h_func_I(xk(8), -1*epsilon, [x_hat_LOS(1), x_hat_LOS(2), xk(1)]', ...
                        [sv_pos_LOS(1), sv_pos_LOS(2), sv_pos_LOS(3)]', atan2(Q_E_filt,I_E_filt));
    h_Q_E = h_func_Q(xk(8), -1*epsilon, [x_hat_LOS(1), x_hat_LOS(2), xk(1)]', ...
                        [sv_pos_LOS(1), sv_pos_LOS(2), sv_pos_LOS(3)]', atan2(Q_E_filt,I_E_filt));
    h_I_L = h_func_I(xk(8), epsilon, [x_hat_LOS(1), x_hat_LOS(2), xk(1)]', ...
                        [sv_pos_LOS(1), sv_pos_LOS(2), sv_pos_LOS(3)]', atan2(Q_L_filt,I_L_filt));
    h_Q_L = h_func_Q(xk(8), epsilon, [x_hat_LOS(1), x_hat_LOS(2), xk(1)]', ...
                        [sv_pos_LOS(1), sv_pos_LOS(2), sv_pos_LOS(3)]', atan2(Q_L_filt,I_L_filt));

    % concatenate into observation matrix 
    h_star_tmp = [h_I_E;h_Q_E;h_I_L;h_Q_L];
    h_star_tmp = Dy^-1*h_star_tmp*Dx; % scale jacobian 
    hx_star = h_star_tmp * [xk(1), xk(3), xk(5), xk(6), xk(8)]';
    % measurement does not consist of v_r, dtrdot, or b
    % add zeros to accomodate this
    h_star = [h_star_tmp(:,1), zeros(4,1), h_star_tmp(:,2), zeros(4,1), ...
                h_star_tmp(:,3:4), zeros(4,1), h_star_tmp(:,5)];    
    % get nominal measurement
    % gen_nominal_signal = @(a, epsilon, dtheta, isCosine)
    ynom_I_E = gen_nominal_signal(a, -1*epsilon, atan2(Q_E_filt,I_E_filt), 1);
    ynom_Q_E = gen_nominal_signal(a, -1*epsilon, atan2(Q_E_filt,I_E_filt), 0);
    ynom_I_L = gen_nominal_signal(a, epsilon, atan2(Q_L_filt,I_L_filt), 1);
    ynom_Q_L = gen_nominal_signal(a, epsilon, atan2(Q_L_filt,I_L_filt), 0);

    % concatenate into nominal measurement 
    y_nom    = [ynom_I_E;ynom_Q_E;ynom_I_L;ynom_Q_L];

    % form adjusted measurement 
    y_star = Dy^-1*(y_filt - y_nom) + hx_star; % jacobian scaling changes measurement 

    % run measurement update
    K = Pk * h_star' * inv(h_star * Pk * h_star' + R);
    innov = y_star - h_star*xk;
    xhatk = xk + K*(innov);
    P_ch = (eye(size(Pk)) - K*h_star);
    Phatk = P_ch*Pk;
    
    % time update 
    % take state vector r_r, v_r, dtr, dtrdot, I, T, b, P
    % control inputs are r_dot_nominal, v_nominal, v_dot_nominal, and gadj
    xk = A*xhatk + B*[xk(2);xk(2);0;g_adj];
    Pk = A*Phatk*A' + Q;
    Pk_r = [Pk_r, sqrt(abs(Pk(1,1)))];
    x_hat_LOS(3) = xk(1);
    range = norm(sv_pos_LOS-x_hat_LOS); 
    range_l = [range_l, range];

    % get next data block 
    [rawSignal, samplesRead] = fread(fid, ...
    dataAdaptCoeff*blksize, settings.dataType);
    % grab I and Q components 
    if mod(i,20) == 0
       counter = counter + 1;
       nb_now = nb(counter);  
    end 
    rawSignal = nb_now*rawSignal';
    rawSignal1=rawSignal(1:2:end);
    rawSignal2=rawSignal(2:2:end);
    rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector

    % calc new time of travel 
    % tau = 1/c*(range + xk(5) + xk(6) + c*(dtrk-dtsv))*1000;
    tau = 1/c*(range + xk(5) + xk(6) + c*(xk(3)-dtsv))*1000;
    rem_tau = tau - floor(tau);

end

data_plt = Pk_r;
kernel = 1/2*ones(1,2);
smoothedData = conv(data_plt, kernel, 'same');


plot(smoothedData(1:end-1))
xlabel('time [ms]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\hat{\sigma_{r_r}}$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
title('Range Covariance vs. time', 'Interpreter', 'latex');
grid on; 
set(gca, 'FontSize', 14);
% Adding a textbox in the upper right corner
annotation('textbox', [0.5, 0.8, 0.2, 0.1], 'String', 'C/N0 ranges from -4 to 4 dB Hz', ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white');

% Save as a PNG with 300 dots per inch (DPI)
print('myHighResPlot', '-dpng', '-r300');

kenerl3= 1/10*ones(1,10);
cn0_list(1000:end) = 2*conv(cn0_list(1000:end), kenerl3,'same'); 
