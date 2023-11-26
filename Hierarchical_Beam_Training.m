% % Hierarchical Beam Training
clc; clear; close all;

%% Transmitter and Receiver
N_Rx = 8;
d_Rx = 0.5;

N_Tx = 16;
d_Tx = 0.5; % half wavelength


%% Receiver with DFT codebook
if mod(N_Rx, 2) == 0
    u_Rx = -1 : 2/N_Rx : 1-(2/N_Rx);
else
    u_Rx = -1 + 1/N_Rx : 2/N_Rx : 1;
end

w_Rx = zeros(N_Rx, N_Rx);
% RF codebook
for k=1:length(u_Rx)
    w_Rx(:,k) = sqrt(1/N_Rx) * exp(-1i*2*pi*d_Rx*(0:N_Rx-1)*u_Rx(k));
end

%% Transmitter with multi-level codebook
if mod(N_Rx, 2) == 0
    u_Rx = -1 : 2/N_Rx : 1-(2/N_Rx);
else
    u_Rx = -1 + 1/N_Rx : 2/N_Rx : 1;
end

w_Rx = zeros(N_Rx, N_Rx);
% RF codebook
for k=1:length(u_Rx)
    w_Rx(:,k) = sqrt(1/N_Rx) * exp(-1i*2*pi*d_Rx*(0:N_Rx-1)*u_Rx(k));
end

%% Channel
AoA = 30;
AoD = 25;
alpha = 1; % channel gain;

% Array Steering Vector
a_Rx = sqrt(1/N_Rx) * exp(-1i*2*pi*d_Rx*(0:N_Rx-1)*sind(AoA)).';
% a_Tx = sqrt(1/N_Tx) * exp(-1i*2*pi*d_Tx*(0:N_Tx-1)*sind(AoD)).';

% number of levels, every new level has 4 times the elements
n = zeros(floor(log(N_Tx)/log(4)), 1); % log4(N_Tx)
u_Tx = zeros(N_Tx, length(n));
a_Tx = zeros(N_Tx, length(n));
H = zeros(N_Rx, N_Tx, length(n));
w_Tx = zeros(N_Tx, N_Tx, length(n));


for i = 1:length(n)
    n(i) = 4^i;
    u_Tx(1:n(i), i) = -1 : 2/n(i) :1-(2/n(i)); % steps
    for k = 1:n(i)
        w_Tx(1:n(i), k, i) = sqrt(1/n(i))*exp(-1i*2*pi*d_Tx*(0:n(i)-1)*u_Tx(k,i)); % Tx beamformer
    end
end


c = zeros(length(u_Rx), N_Tx, length(n));

for i = 1:length(n)
    a_Tx(1:n(i), i) = sqrt(1/n(i))*exp(-1i*2*pi*d_Tx*(0:n(i)-1)*sind(AoD)).';
    H(:, 1:n(i), i) = sqrt(1/(n(i)*N_Rx))*alpha*a_Rx*a_Tx(1:n(i),i)';
    if i==1
        for j=1:n(i)
            for k=1:length(u_Rx)
                c(k,j,i) = w_Rx(:,k)'*H(:,1:n(i),i)*w_Tx(1:n(i),j,i);
            end
        end
    else
        p = 1 + (I_col(i-1)-1)*4;
        if p <= 2
            p_min = 1;
        else
            p_min = p-2;
        end
            
        for j = p_min : p+2
            for k = 1:length(u_Rx)
                c(k,j,i) = w_Rx(:,k)'*H(:,1:n(i),i)*w_Tx(1:n(i),j,i);
            end
        end
    end

    temp = c(:, 1:n(i), i);
    Pow = abs(temp);
    
    %[M(i), Idx(i)] = max(Pow(:));
    %[I_row(i), I_col(i)] = ind2sub(size(Pow), Idx(i)); 
    [I_row(i), I_col(i)] = find(Pow == max(max(Pow)));
    
    AoD_est = asind(u_Tx(I_col(i), i));
    
    figure();
    bar3(Pow); % 'detached'
    
    str=sprintf('Level (%d) Beam Training', i);
    title(str);
    xlim([0, n(i)+1]);
    ylim([0, N_Rx+1]);
    xlabel('Transmit Beams');
    ylabel('Receive Beams');
    
end


