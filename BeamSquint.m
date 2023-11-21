clear;
clc;
close all;
%% Basic Electromagnetic Parameters
Frequency = 28e9;
Lightspeed = physconst('LightSpeed');

%% Array Parameters
N = 16;
% X = (1:N)*WaveLength/2;
% alpha = zeros(1,N);

%% ArrayFactor Samping
Ns =500;% Sampling number
theta = linspace(-90,90,Ns);
E =zeros(3,Ns);

doa_theta = -50;

% mutli-band
freq = [27, 28, 29]*1e9;
% Wavelength = Lightspeed./freq;
% Wavenumber = 2*pi./Wavelength;

% spacing = (1:N)*Wavelength(1)/2;
W = exp(-1i*((0:N-1)'*pi*sind(doa_theta)));

for k=1:length(freq)
    for num = 1:Ns
        E(k, num)=W'*exp(-1i*((0:N-1)'*pi*sind(theta(num))*(freq(k)/Frequency)));
    end
end


%% plot figure
figure(1);
plot(theta, db(E(1,:)),'LineWidth',2);
hold on;
plot(theta, db(E(2,:)),'LineWidth',2);
plot(theta, db(E(3,:)),'LineWidth',2);
grid on;
xline(doa_theta,'m--','LineWidth',2);
xlim([-90, 90]);
% ylim([-100, 50]);
xticks([-90:30:90]);
xlabel('\theta(\circ)');
ylabel('dB');
set(gca,'Fontsize',13)
legend('27 GHz', '28 GHz', '29 GHz');

% figure(2);
% plot(theta, db(E(1,:))-max(db(E(1,:))),'LineWidth',2);%normalized
% grid on;
% hold on;
% plot(theta, db(E(2,:))-max(db(E(2,:))),'LineWidth',2);
% plot(theta, db(E(3,:))-max(db(E(3,:))),'LineWidth',2);
% xlim([-90, 90]);
% % ylim([-100, 0]);
% xticks([-90:30:90]);
% xlabel('\theta(\circ)');
% ylabel('dB');
% set(gca,'Fontsize',13)
% legend('27 GHz', '28 GHz', '29 GHz');