%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% SIMPLE UNIFORM LINEAR ARRRAY 
% WITH VARIABLE NUMBER OF ELEMENTS 
% MATRIX IMPLEMENTATION 
% COPYRIGHT RAYMAPS (C) 2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

f = 28e9;     % carrier frequency
c = 3e8;      % light speed
lambda = c/f; % wavelength
d = lambda/2; 
no_elements = 8; 
theta = 0 : pi/180 : 2*pi;

n = 1:no_elements; 
n = transpose(n); 
A = (n-1)*(1i*2*pi*d*cos(theta)/lambda); %sin(theta)
X = exp(-A); 
w = ones(1, no_elements); 
r = w*X;
r_dB = abs(r);

figure(1);
polarplot(theta, r_dB, 'b', 'LineWidth', 2);
thetalim([0, 180]); 
hold on;

ref = max(r_dB)*ones(size(theta));
polarplot(theta, ref, 'r-.', 'LineWidth', 2);

title ('Gain of a Uniform Linear Array');
set(gca,'fontsize', 12);

 
