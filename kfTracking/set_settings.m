% settings.m
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

% get jacobians before tracking starts
syms phi theta psi % euler angles for jacobian computation
Qbe_inv_sym = compute_Qbe_inv(phi, theta, psi);
R_BN_sym = enuToBody(phi, theta, psi);

tensor = sym(zeros(3,3,3));
tensor2 = sym(zeros(3,3,3));
for i = 1:3
    for j = 1:3
        for k = 1:3
            tensor(i,j,k) = diff(Qbe_inv_sym(i,j), E(k));
            tensor2(i,j,k) = diff(R_BN_sym(i,j), E(k));
        end
    end
end

dQbe_inv_dE = reshape(tensor, 9, 3);
dR_BN_dE = reshape(tensor2, 9, 3);

% Define anonymous functions for dQ_be and dR_BN
settings.dQ_be = @(phi_val, theta_val, psi_val) double(subs(dQbe_inv_dE, [phi, theta, psi], [phi_val, theta_val, psi_val]));
settings.dR_BN = @(phi_val, theta_val, psi_val) double(subs(dR_BN_dE, [phi, theta, psi], [phi_val, theta_val, psi_val]));
%dQ_be_now = settings.dQ_be(1, 2, 3); % Replace 1, 2, 3 with your actual values for phi, theta, and psi

