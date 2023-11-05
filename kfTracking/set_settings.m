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

% fetch continuous time IMU model 
[Ac, Bc, Gc] = getIMUModelContinuous();

% Discretize the system
settings.TsIMU = 0.01;  % Sample time (adjust as needed)

% c2d is not defined for syms, so I'll do it myself
Ad = expm(Ac*settings.TsIMU);

if(det(A) == 0)
    syms tau 
    Bd = vpa(int(expm(Ac*tau)*Bc,tau,[0 settings.TsIMU]));
else
    Bd = Ac\(Ad-eye(size(Ad)))*Bc;
end

% add remaining 2 state clock, iono, tropo, and signal amplitude models 
% dtr_{k} = dtr_{k-1} + Ts*dtr_dot_{k-1} + 0.5*normrnd(0, sig_d)*Ts^2
% dtr_dot_{k} = dtr_dot_{k-1} + normrnd(0, sig_d)
% I_{k} = I_{k-1} + (-1/tau_I * I_{k-1} + normrnd(0, sig_I))*Ts
% T_{k} = T_{k-1} + (-1/tau_T * T_{k-1} + normrnd(0, sig_T))*Ts
% a = a + normrnd(0, sig_a)*Ts

% new state vector (20 total states)
%[r v E ba bg dtr dtr_dot I T a}]

% allocate space for non-IMU models 
[an,am] = size(Ad); % 9 state PVA (should be 9x15)
[bn,bm] = size(Bd); 
[gn,gm] = size(Gc);
full_Ad = zeros(an+6+5,am+5); % 9 state PVA, 6 biases, and 5 above models 
full_Bd = zeros(bn+6+5,bm);
full_Gc = eye(gn+5+6, gm+5+6);
full_Ad(1:an,1:am) = Ad;
full_bd(1:bn,:) = Bd;
full_Gc(1:gn,1:gm) = Gc; 

% Clock TXCO
sig_d = 5e-10^2; % Discrete clock drift noise variance

% amplitude 
sig_a = 5; 

% Iono and tropo stochastic parameters 
tau_I = 30*60;    % About a 30-minute time constant
tau_T = 2*3600;   % About a 2-hour time constant
sig_I = 2^2;        % Should be about between 1-5 meters
sig_T = .5*sig_I^2; % Should be roughly 1/2 of Iono

% process noise information for IMU, iono, tropo, signal amplitude
% TODO set into calls like "tactical, navigation, auto to fetch this 
Q = diag([0 0 0 57e-6 57e-6 57e-6 57e-6 57e-6 57e-6 1e6 1e6 1e6 1e6 1e6 ...
            1e6 sig_d*.5 sig_d sig_I sig_T sig_a]);

Qd = full_Gc*Q*full_Gc*settings.TsIMU; 

% put back in IMU
full_Qd = Qd;

% FOGMP biases accelerometer/gyro 
tau_b_acc = 3600;
tau_b_gy = 3600;
A_biases = [(1-settings.TsIMU/tau_b_acc)*eye(3) zeros(3,3); 
             zeros(3,3) (1-settings.TsIMU/tau_b_gy)*eye(3)];
full_Ad(an+1:an+1+5,an+1:an+1+5) = A_biases;

% 2 state clock model 
A_2sc = [1 settings.TsIMU; 0 1];
full_Ad(an+6+1:an+6+2,an+6+1:an+6+2) = A_2sc;

% iono and tropo
A_it = [(1-settings.TsIMU/tau_I) 0; 0 (1-settings.TsIMU/tau_T)];
full_Ad(an+6+3:an+6+4,an+6+3:an+6+4) = A_it;

% amplitude 
full_Ad(end,end) = 1;

% syms nu
% Gd = vpa(int(expm(A*tau)*G,tau,[0 Ts]));
% big=expm([A G*q*G;zeros(3) -A']*Ts);
% Ad=big(1:3,1:3);
% Qd=big(1:3,4:6) * Ad';