% % Exhaustive Beam Training
clc; clear; close all;

%% Transmitter
N_Tx = 8;
d_Tx = 0.5; % half wavelength

if mod(N_Tx, 2) == 0
    u_Tx = -1 : 2/N_Tx : 1-(2/N_Tx);
else
    u_Tx = -1 + 1/N_Tx : 2/N_Tx : 1;
end

% RF codebook
for k=1:length(u_Tx)
    w_Tx(:,k) = sqrt(1/N_Tx) * exp(-1i*2*pi*d_Tx*(0:N_Tx-1)*u_Tx(k));
end


%% Receiver
N_Rx = 8;
d_Rx = 0.5;

if mod(N_Rx, 2) == 0
    u_Rx = -1 : 2/N_Rx : 1-(2/N_Rx);
else
    u_Rx = -1 + 1/N_Rx : 2/N_Rx : 1;
end

% RF codebook
for k=1:length(u_Rx)
    w_Rx(:,k) = sqrt(1/N_Rx) * exp(-1i*2*pi*d_Rx*(0:N_Rx-1)*u_Rx(k));
end


% Angles
AoA = 30;
AoD = 57;

alpha = 1; % channel gain;

% Array Steering Vector
a_Tx = sqrt(1/N_Tx) * exp(-1i*2*pi*d_Tx*(0:N_Tx-1)*sind(AoD)).';
a_Rx = sqrt(1/N_Rx) * exp(-1i*2*pi*d_Rx*(0:N_Rx-1)*sind(AoA)).';

% Channel
H = sqrt(N_Tx*N_Rx)*alpha*a_Rx*a_Tx';

%% Beam Training
c = zeros(N_Rx, N_Tx);
for j = 1:length(u_Tx)
    for k = 1:length(u_Rx)
        c(k,j) = w_Rx(:,k)'*H*w_Tx(:,j);
    end
end

C = abs(c);
C_vec = C(:); % vectorize
[Entry, Idx] = max(abs(C_vec));
[I_row, I_col] = ind2sub(size(C), Idx);


Tx_beam_number = I_col; % matches with AoD
Rx_beam_number = I_row; % matches with AoA

figure(1);
% set(0, 'DefaultAxesFontSize', 14, 'DefaultTextFontSize', 14);
bar3(C); % 'detached'
title('Ideal Beam Training');
xlabel('Transmit Beams');
ylabel('Receive Beams');

figure(2)
stem3 (C, 'Marker', 's', ...
'MarkerEdgeColor' , 'm', ...
'MarkerFaceColor' , 'g') ;
title('Ideal Beam Training');
xlabel('Transmit Beams');
ylabel('Receive Beams');
xlim([1, N_Tx]);
ylim([1, N_Rx]);
hold on;
grid on;

