% % UPA array steering vector
clear; clc;

%% Antenna paras
Nx = 8;
Ny = 8;
N = Nx * Ny;
theta = 40;   % 方向角
phi = 130;     % 俯仰角

nx = 0:Nx-1;
ny = 0:Ny-1;
%% UPA steering vec
% Low-Complexity Separable Beamformers for Massive Antenna Array Systems
% array_y = sqrt(1/Ny)*exp(-1i * pi * ny * sind(theta) * sind(phi)).';
% array_z = sqrt(1/Nx)*exp(-1i * pi * nx * cosd(phi)).';
% Array = kron(array_y, array_z);

Array = zeros(Nx, Ny);
for k = 0 : Nx-1
    for j = 0 : Ny-1
        Array(k + 1, j + 1) = exp(-1i * pi * (k*cosd(phi) + j*sind(theta) * sind(phi)))/sqrt(Nx*Ny);
    end
end



%% 2D DFT matrix/codebook for UPA
DFT_H = zeros(N, N);
% for k = 0 : N-1
%     v = floor(k / Nx);
%     u = k - v * Nx;
%     for j = 0 : N-1
%         y = floor(j / Nx);
%         x = j - y * Nx;
%         DFT_H(k + 1, j + 1) = exp(-1i * 2 * pi * (u * x/ Nx + v * y / Ny));
%     end
% end

% Fx = dftmtx(Nx);
% Fy = dftmtx(Ny);
gy = 0 : -2/Ny : -2*(Ny-1)/Ny;
gy((gy<-1)) = gy((gy<-1)) + 2;
gx = 0 : -2/Nx : -2*(Nx-1)/Nx;
gx((gx<-1)) = gx((gx<-1)) + 2;


Fy = exp(1i * pi * ny' * gy);
Fx = exp(1i * pi * nx' * gx);

DFT_H = kron(Fy.', Fx);

%% Beam training
received_seq = abs(Fx*Array*Fy);
[I_row, I_col] = find(received_seq == max(max(received_seq)));

% https://blog.csdn.net/qq_23947237/article/details/89925088
% [ -(N-1)/2N,  (N-1)/2N ]
dft_phi = gx(I_row);
est_phi = acosd(dft_phi);
dft_theta = gy(I_col)/sind(est_phi);
est_theta = asind(dft_theta);


disp(['estimated phi:', num2str(est_phi)]);
disp(['estimated theta:', num2str(est_theta)]);

figure(1);
bar3(received_seq); % 'detached'




