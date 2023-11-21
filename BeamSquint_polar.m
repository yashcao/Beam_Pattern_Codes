%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% SIMPLE UNIFORM LINEAR ARRRAY 
% WITH VARIABLE NUMBER OF ELEMENTS 
% MATRIX IMPLEMENTATION 
% COPYRIGHT RAYMAPS (C) 2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

f0 = 28e9;     % carrier frequency
c = 3e8;      % light speed
% lambda = c/f; % wavelength
% d = lambda/2; 
no_elements = 16;
Ns = 500;
% theta = 0 : 2*pi/Ns : 2*pi;
theta = linspace(0,2*pi,Ns);

n = 1:no_elements; 
N = transpose(n);
%w = ones(1, no_elements);
theta0 = 50;
W = (n-1)*(1i*pi*cosd(theta0)); %sin(theta)

f = [27, 28, 29]*1e9;
r_dB =zeros(3, Ns);
for k=1:length(f)
    A_scan = (N-1)*(1i*pi*cos(theta)*(f(k)/f0)); %sin(theta)
    A = exp(-A_scan);
    
    r = exp(W)*A;
    r_dB(k,:) = abs(r);
end



figure(1);
polarplot(theta, r_dB(1,:), 'LineWidth', 1.5);
hold on;
polarplot(theta, r_dB(2,:), 'LineWidth', 1.5);
polarplot(theta, r_dB(3,:), 'LineWidth', 1.5);
thetalim([0, 180]);

% polarplot([0, 0], [20, max(max(r_dB))],'-','LineWidth',2,'Color','b')
polarplot([0 0; 0 50]*pi/180, [0 0; 0 1]*max(max(r_dB)), '-.','LineWidth',1.5,'Color','g');
% ref = max(r_dB)*ones(size(theta));
% polarplot(theta, ref, 'r-.', 'LineWidth', 2);

title ('Gain of a Uniform Linear Array');
set(gca,'fontsize', 12);
legend('27 GHz', '28 GHz', '29 GHz', 'Location', 'Best');


 
