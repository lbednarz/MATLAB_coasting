%% Constants, data inputs 

% This script initializes settings for the tracking loop.
addpath Inertial\

% constants
% initials for position and signal quality
settings.c = 2.99792458e8;             % speed of light [m/s]
settings.carrier_freq = 1575.42e6;     % carrier frequency [Hz]
settings.lam = settings.c/settings.carrier_freq; % carrier wavelength [m]
settings.codeFreqBasis = 1/1.023e6;           % chip duration [s]
settings.Ts = 10e6/2;                         % Hz sample rate
settings.Tacc = 1e-3;                         % accumulation period [s] 
settings.f_IF = 0e6;                          % intermediate frequency [Hz]

% data to use 
settings.RFdata = fullfile('..','..','GPSStatic_BW8_gain65_10Msps__REAL.bin');
settings.trackingData = "data\trackingResults96.mat";
settings.priors = "data\extras97.mat";
dataAdaptCoeff = 2; % 1 means real data at 8bit size, 2 is IQ complex data

%% Kalman filter parameters 

syms phi theta psi_var phi_dot theta_dot psi_var_dot fb1 fb2 fb3 lat lon;
% fetch continuous time IMU model
tau_b_acc = 3600;
tau_b_gy = 3600;

% iono and tropo information
tau_I = 30*60;    % About a 30-minute time constant
tau_T = 2*3600;   % About a 2-hour time constant
sig_I = 2;        % Should be about between 1-5 meters
sig_T = .5*sig_I; % Should be roughly 1/2 of Iono

% Discretize the system
settings.TsIMU = 0.01;  % Sample time (adjust as needed)

% Clock TXCO
sig_d = 5e-10; % Discrete clock drift noise variance

% amplitude uncertainty  
sig_a = 5; 

% accelerometer and gyro noise 
sig_acc = 57e-6;
sig_gy = 57e-6;
bias_acc = 1e6;
bias_gy = 1e6;

% process noise information for IMU, iono, tropo, signal amplitude
% TODO set into calls like "tactical, navigation, auto to fetch this 
settings.W = diag([0 0 0 sig_acc*ones(1,3) sig_gy*ones(1,3) ... 
                    bias_acc*ones(1,3) bias_gy*ones(1,3) ... 
                    sig_d*.5 sig_d sig_I sig_T sig_a]).^2;

% Qd = full_Gc*Q*full_Gc*settings.TsIMU; needs to be at runtime  



