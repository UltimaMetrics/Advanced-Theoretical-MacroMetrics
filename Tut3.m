% ------------------------------------------------------------------------------------
% This code implements the Homoskedastic linear Gaussian state space model estimated 
% the Kalman Filter and Smoother using the Durbin and Koopman (2002) algorithm for 
% state-space models.
% ************************************************************************************
% The local level model is:
%
%     Y(t) = B(t) + u(t) 
% 
%  with u(t)~N(0,H).
% The state equation is
%
%            B(t) = B(t-1) + error
%
%

clear all;
clc;
%----------------------------------LOAD DATA----------------------------------------
% Load data from the second lab session
load lab2.dat -ascii
t = lab2(:,1);
y = lab2(:,2:4);
% Recall that y containts unemployment, interest rates and inflation
T = size(y,1);

% Create a vector of ones
ii = ones(T,1);
% plot(y);
% Demean and standardize data
stdy = std(y);
ys = (y - mean(y,1))./stdy;
y1 = y(:,1);
s2 = 0.85;%e
se2 = 0.005;  %n
[bhatll,Vtt] = KalFilt(y1,ii,ii*s2,ii*se2,1,1,T,1,1);
' Local level model estimated'

% Now plot the data with the local level model
plot(t,[y1 bhatll'])
