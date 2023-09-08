% load in prior infomation from PLL/DLL
% for this demo, we are interested in channels 2,3,8,9
close all; clear; clc; 
addpath 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GNSS_SDR-master'
addpath data

load("data\extras97.mat"); % should set a variable called "extra_out"
load("data\trackingResults96.mat")
dataAdaptCoeff = 1; 
channel_to_track = 8; % channel we want to track
PRN = 13; % PRN in channel
frame_starts = extra_out.sample_num; % contains byte start of each
                                     % subframe - 5515
% should be 220614790
byte_frame_start = trackResults(channel_to_track).absoluteSample(frame_starts(1,channel_to_track));
settings = initSettings();

% get RF file to track 
fileName = 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GPSStatic_BW8_gain65_10Msps__REAL.bin';
[fid, message] = fopen(fileName, 'rb');

% Move to location of first prior 
  fseek(fid, ...
    dataAdaptCoeff*byte_frame_start, ...
    'bof');

remCodePhase = 0;
codeFreq      = settings.codeFreqBasis;
codePhaseStep = codeFreq / settings.samplingFreq;
blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
[rawSignal, samplesRead] = fread(fid, ...
    dataAdaptCoeff*blksize, settings.dataType);

rawSignal = rawSignal';

% get navbits 
nb = extra_out.extra_nb{1,channel_to_track};
nb(nb == 0) = -1; % encode NB into 1s and -1s 

% raw signal w/o navbit
rawSignal = nb(1)*rawSignal;

% repeat for each measurement epoch. Next batch will get a new navbit that
% will be applicable for 20 code periods.

%% build priors 
% initials for position and signal quality
c = 2.99792458e8;         % speed of light in m/s
L1 = 1575.42e6;           % L1 carrier frequency
lam = c/L1;               % carrier wavelength
Tc = 1/1.023e6;           % chip duration

priors_1 = extra_out.extras_epoch.epoch1;
% intial estimates 
%[x_hat, y_hat, z_hat, cdtr_hat] = priors_1.state; 
x_hat = priors_1.state;
range = norm(priors_1.satpos(:,channel_to_track)-x_hat(1:3));
tau = priors_1.tot*1000; % travel time in ms 
rem_tau = tau - floor(tau); % retrieve just the time into current period 
epsilon = .5*Tc;
dtr = x_hat(4)/c;
dtsv = priors_1.dtsv(8);
Tr = priors_1.tropo; % tropo from SDR is high and reports no covariance 
Io = 0; % because T is so high, we will assert 0 tropo with some covariance
R_EB = LOSRotationMatrix(x_hat(1:3), priors_1.satpos(:,channel_to_track));
% initial covairances
P0 = priors_1.Q; % assumes states are [x, y, z, c*dtr]
P0(1:3,1:3) = R_EB*P0(1:3,1:3); % rotate into body frame
% conservatively grab 1D pos and clock uncertainty 
P0 = [P0(2,2), P0(2,4); P0(4,2), P0(4,4)].*10;
% fill out P0 to include tropo and iono



%% measurement models 
% get spreading code for specified sample rate, PRN, and time of travel
% using get_x(rem_tot, PRN, Ts)
% x comes from get_x function
% isCosine allows for making I and Q replicas from 1 function 
gen_local_replica = @(x, t, lambda, range, I, T, c, dtr, dtsv, isCosine) ...
    x .* (isCosine .* cos(2*pi/lambda * (range - I + T + c*(dtr - dtsv))) + ...
    (1-isCosine) .* sin(2*pi/lambda * (range - I + T + c*(dtr - dtsv))));

% Define R(epsilon) 
R_func = @(epsilon) (abs(epsilon) < Tc) .* (1 - abs(epsilon)/Tc);

% Define R'(epsilon)
R_prime_func = @(epsilon) (epsilon < 0) .* 1.023e6 + ...
                          (epsilon > 0) .* -1.023e6;

r_los = @(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3) ...
    ([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]) / ...
    norm([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]);

h_func_I = @(a, epsilon, r_r1, r_r2, r_r3, v_r1, v_r2, v_r3, r_sv1, r_sv2, r_sv3, v_sv1, v_sv2, v_sv3) ...
    [R_func(epsilon), ... 
    a/c * R_prime_func(epsilon) * dot(r_los(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3), [v_r1; v_r2; v_r3] - [v_sv1; v_sv2; v_sv3]), ...
    a * R_prime_func(epsilon) / c, ... 
    a * R_prime_func(epsilon) / c, ... 
    a * R_prime_func(epsilon) / c];

h_func_Q = @(a, epsilon, r_r1, r_r2, r_r3, v_r1, v_r2, v_r3, r_sv1, r_sv2, r_sv3, v_sv1, v_sv2, v_sv3) ...
    [0, ... 
    2*pi*a/lam * R_func(epsilon) * dot(r_los(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3), [v_r1; v_r2; v_r3] - [v_sv1; v_sv2; v_sv3]), ...
    -2*pi*a/lam * R_func(epsilon), ... 
    2*pi*a/lam * R_func(epsilon), ... 
    2*pi*a/lam * R_func(epsilon)];


%% dynamic models 

% Initialization
sim_time = 1;
Ts = (10*1.023e6)^-1; % Hz sample rate

% Initial values for position and signal quality
c = 2.99792458e8; % Speed of light in m/s
L1 = 1575.42e6; % L1 carrier frequency
lam = c/L1; % Carrier wavelength
CN0dBHz = -10; % Carrier to noise density ratio dBHz
CN0 = 10^(CN0dBHz/10); % Carrier to noise density ratio
R = 1/CN0; % Variance of measurement noise
sig_y = sqrt(R); % SD of measurement noise
PdBHz = -130; % Signal power dB/Hz
P = 10^(PdBHz/10);
n = 1; % Initial noise
r_r = 0; % On earth's surface ENU
v_r = 1; % m/s
r_sv = 22000e3; % Meters ENU
dt_sv = 1e-3;
sig_v = .1; % m/s

% Iono and tropo
tau_I = 30*60; % About a 30-minute time constant
tau_T = 2*3600; % About a 2-hour time constant
sig_I = 2; % Should be about between 1-5 meters
sig_T = .5*sig_I; % Should be roughly 1/2 of Iono
I = 1; % Meters
T = .5; % Meters

% Accelerometer with bias
sig_a = 57e-6; % SD accelerometer noise m/s^1.5
sig_b = 14e-6; % Continuous SD of accelerometer bias RW noise
tau_b = 3600; % Time constant of accelerometer bias RW
b = 14e-6; % Initial guess at bias

% Clock TXCO
sig_d = 5e-8 * sqrt(Ts); % Discrete clock drift noise variance
b_c = 1e-8; % Initial receiver clock bias
d_c = 2*sig_d; % Initial receiver clock drift

% take state vector r_r, v_r, dtr, dtrdot, I, T, b 
% r_eq = r_{k-1} + (v_eq + (r_dot_nominal_{k-1} - v_nominal_{k-1}))*Ts
% v_eq = v_k = v_{k-1} + (-1*b + (v_dot_nominal_{k-1}-g_adj) - normrnd(0, sig_a))*Ts
% dtr_{k} = dtr_{k-1} + Ts*dtr_dot_{k-1} + 0.5*normrnd(0, sig_d)*Ts^2
% dtr_dot_{k} = dtr_dot_{k-1} + normrnd(0, sig_d)
% I_{k} = I_{k-1} + (-1/tau_I * I + normrnd(0, sig_I))*Ts
% T_{k} = T_{k-1} + (-1/tau_T * T + normrnd(0, sig_T))*Ts
% b_{k} = b_{k-1} + (-1/tau_b * b_{k-1} + normrnd(0, sig_b))*Ts

% Assuming values for Ts, tau_I, tau_T, tau_b, and g_adj are provided:

A = [1, Ts, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, -Ts;
     0, 0, 1, Ts, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 1 - Ts/tau_I, 0, 0;
     0, 0, 0, 0, 0, 1 - Ts/tau_T, 0;
     0, 0, 0, 0, 0, 0, 1 - Ts/tau_b];

B = [Ts, -Ts, 0;
     -Ts, Ts, -Ts;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0];

Q = diag([0, sig_a^2, 0.5^2 * sig_d^2 * Ts^2, sig_d^2, sig_I^2, sig_T^2, sig_b^2]);
