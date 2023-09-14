clear; close all; clc;
F = eye(8);
G = eye(8);
H = eye(8);
Qsimple = eye(8);
Rsimple = 1*eye(8);
sysSimple = ss(F, G, H, []);
[~,~,Psimple] = kalman(sysSimple, Qsimple, Rsimple);
