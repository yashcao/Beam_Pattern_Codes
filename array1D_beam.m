clear;
clc;
close all;
%% Basic Electromagnetic Parameters
Frequency = 10e9;
Lightspeed = physconst('LightSpeed');
Wavelength = Lightspeed/Frequency;
Wavenumber = 2*pi/Wavelength;

%% Array Parameters
N = 16;
X = (1:N)*Wavelength/2;
W =  ones(1,N);
alpha = zeros(1,N);

%% ArrayFactor Samping
Ns =500;% Sampling number
theta = linspace(-90,90,Ns);
E =zeros(1,Ns);

W = exp(1j*(Wavenumber*X*sind(-50)+alpha));

for num = 1:Ns
    E(num)=W*exp(1j*(Wavenumber*X*sind(theta(num))+alpha))';
end
%% plot figure
figure(1);
plot(theta,db(E),'LineWidth',2);
grid on;
xlim([-90, 90]);
xticks([-90:30:90]);
xlabel('\theta(\circ)');
ylabel('dB');
set(gca,'Fontsize',13)

figure(2);
plot(theta, db(E)-max(db(E)),'LineWidth',2);%normalized
grid on;
xlim([-90, 90]);
xticks([-90:30:90]);
xlabel('\theta(\circ)');
ylabel('dB');
set(gca,'Fontsize',13)
