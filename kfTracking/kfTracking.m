close all; clear; clc; 

% load in prior infomation from PLL/DLL
% for this demo, we are interested in channels 2,3,8,9
addpath 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GNSS_SDR-master'
load("extras96.mat"); % should set a variable called "extra_out"
dataAdaptCoeff = 1; 
channel_to_track = 8;
frame_starts = extra_out.sample_num; % contains byte start of each
                                     % subframe 
settings = initSettings();

% get RF file to track 
fileName = 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\GPSStatic_BW8_gain65_10Msps__REAL.bin';
[fid, message] = fopen(fileName, 'rb');

% Move to location of first prior 
  fseek(fid, ...
    dataAdaptCoeff*frame_starts(1,channel_to_track), ...
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
priors_1 = extra_out.extras_epoch.epoch1;
%[x_hat, y_hat, z_hat, cdtr_hat] = priors_1.state; 
x_hat = priors_1.state;

%% measurement models 

% initials for position and signal quality
c = 2.99792458e8;         % speed of light in m/s
L1 = 1575.42e6;           % L1 carrier frequency
lam = c/L1;               % carrier wavelength
Tc = 1/1.023e6;           % chip duration
test_epsilon = 0.5 * Tc;  % For example

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
sig_I = 1.5; % Should be about between 1-5 meters
sig_T = .5; % Should be roughly 1/2 of Iono
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

state_out = [];

for i = 1:sim_time
    % Propagate states
    I = I + Ts*(-1/tau_I * I + normrnd(0, sig_I)); % Ionospheric position offset
    T = T + Ts*(-1/tau_T * T + normrnd(0, sig_T)); % Tropospheric position offset
    b = b + Ts*(-1/tau_b * b + normrnd(0, sig_b)); % Accelerometer bias
    v_r = v_r + normrnd(0, sig_v); % Velocity
    r_r = r_r + Ts*v_r; % New position
    w_d = normrnd(0, sig_d); % Clock noise
    b_c = b_c + d_c*Ts + 0.5*w_d*Ts^2; % Clock bias
    d_c = d_c + w_d; % Clock drift
    dr = r_r - r_sv;
    r = norm(dr); % r_r is the 3D user position, r_sv is 3D space vehicle position from ephemeris
    %state_out = [state_out; I, T, b, v_r, r_r];
end

% Define the symbols
% syms a epsilon R R_prime real
% syms r_r [3 1]
% syms v_r [3 1]
% syms r_sv [3 1]
% syms v_sv [3 1]

% shared line of sight vector
% dr = (r_r - r_sv);
% r_los = (r_r - r_sv) / norm(dr);

% % Define the terms for y_I and y_Q
% h_lin_I = [R; a/c * R_prime * dot(r_los, v_r - v_sv); a * R_prime / c; a * R_prime / c; a * R_prime / c];
% h_lin_Q = [0; 2*pi*a/lam * R * dot(r_los, v_r - v_sv); -2*pi*a/lam * R; 2*pi*a/lam * R; 2*pi*a/lam * R];
% 
% % Convert to functions
% variables = [a, R, R_prime, r_r.', v_r.', r_sv.', v_sv.'];
% h_func_I = matlabFunction(h_lin_I, 'Vars', variables);
% h_func_Q = matlabFunction(h_lin_Q, 'Vars', variables);

% r_los = @(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3) ...
%     ([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]) / ...
%     norm([r_r1; r_r2; r_r3] - [r_sv1; r_sv2; r_sv3]);
% 
% h_func_I = @(a, R, R_prime, r_r1, r_r2, r_r3, v_r1, v_r2, v_r3, r_sv1, r_sv2, r_sv3, v_sv1, v_sv2, v_sv3) ...
%     [R; a/c * R_prime * dot(r_los_anon(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3), [v_r1; v_r2; v_r3] - [v_sv1; v_sv2; v_sv3]); ...
%     a * R_prime / c; a * R_prime / c; a * R_prime / c];
% 
% h_func_Q = @(a, R, R_prime, r_r1, r_r2, r_r3, v_r1, v_r2, v_r3, r_sv1, r_sv2, r_sv3, v_sv1, v_sv2, v_sv3) ...
%     [0; 2*pi*a/lam * R * dot(r_los_anon(r_r1, r_r2, r_r3, r_sv1, r_sv2, r_sv3), [v_r1; v_r2; v_r3] - [v_sv1; v_sv2; v_sv3]); ...
%     -2*pi*a/lam * R; 2*pi*a/lam * R; 2*pi*a/lam * R];
% 
% % define R and R' 
% Tc = 1/1.023e6; % Define your Tc value here
% 
% % Define R(epsilon) 
% R_func = @(epsilon) (abs(epsilon) < Tc) .* (1 - abs(epsilon)/Tc);
% 
% % Define R'(epsilon)
% R_prime_func = @(epsilon) (epsilon < 0) .* 1.023e6 + ...
%                           (epsilon > 0) .* -1.023e6;
