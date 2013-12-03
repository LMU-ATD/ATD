% Generation of Sample Radar Signal and a plot of SNR vs Amplitude with varying noise power
%   MATLAB m-file
%   Author: Terry Kong
%   Radar Signal Characteristics:
%       Range(m) = 92600.00, Range(nm) = 50.00
%       Resoluion(dR or deltaR)[in meters] = 154.33
%       Number of Cells = 600
%       Pulse Width(us) = 1.03
%       Interpulse Period(ms) = 0.60
%       Signal Amplitude = 100.00
%       Standard Deviation of Noise = 1.00
%   Variables: 
%       Anbins = number of bins for plots, 
%       rangemax = max range considered,
%       ncells = number of cells in ATD, 
%       Tinter = interpulse period
%       timemax = max time period looked at, 
%       sig_amp = target amplitude(not considering range)
%       sig = standard deviation of [Gaussian] noise (related to power of noise)

clear all; clc; figure(1); clf;
%Constants
nm = 1852; % 1 nautical mile = 1852 m
c = 299279458; % speed of light
mu = 0;
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
nbins = 50;
rangemax = 50*nm; 
ncells = 600;
Tinter = 0.0006; % Interpulse period
timemax = 0.025;
sig_amp = 100; % signal amplitude
sig = 1;
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% delta(R)
dR = rangemax/ncells;
% Pulse width
tau = dR*2/(c); 
% t = time observed
t = 0:tau:timemax; 
% D = delay of pulses
Delay = 0:Tinter:timemax;
% number of points
num = length(t);

% Using Normal Distributions
mu_mat = zeros(1,num) + mu;
sig_mat = zeros(1,num) + sig;
inphase = normrnd(mu_mat,sig_mat);
%A = 100;
%inphase = inphase + A;
quad = normrnd(mu_mat,sig_mat);

% Adding Signal to Noise
signal = zeros(1,length(inphase));
i = 1;
for j = 1:length(t)
    if t(j) > Delay(i)
        signal(j) = signal(j) + sig_amp;
        i = i + 1;
        if i > length(Delay)
            break;
        end
    end
end
signalPower = mean(signal.^2);
noisePower = mean(inphase.^2 + quad.^2);
inphase = inphase + signal;
% plotting signal + noise
subplot(4,1,1)
plot(t,inphase); title(sprintf('Inphase = N_x(0.0,%0.1f) + Sig(w/ amp = %0.2f)',sig^2,sig_amp))
xlabel('t'); ylabel('Amplitude(V)')
subplot(4,1,2)
plot(t,quad); title(sprintf('Quad = N_y(0.0,%0.1f)',sig^2))
xlabel('t'); ylabel('Amplitude(V)')
subplot(4,1,3)
sig_plus_noise = sqrt(inphase.^2 + quad.^2);
plot(t,sig_plus_noise); title(sprintf('Signal plus Noise (Amplitude) [SNR = %0.2fdb]',10*log10(signalPower/noisePower)))
xlabel('t'); ylabel('Amplitude(V)')
subplot(4,1,4)
% bins
bins = linspace(0,max(sig_plus_noise),nbins*1000);
ksdensity(sig_plus_noise,bins);
title('Rician Distribution'); xlim([0,max(sig_plus_noise)])
xlabel('Amplitude'); ylabel('Probability');
fprintf('Range(m) = %0.2f, Range(nm) = %0.2f\n',rangemax, rangemax/nm);
fprintf('Resoluion(dR or deltaR)[in meters] = %0.2f\n',dR);
fprintf('Number of Cells = %d\n',ncells);
fprintf('Pulse Width(us) = %0.2f\n',tau*1000000);
fprintf('Interpulse Period(ms) = %0.2f\n',Tinter*1000);
fprintf('Signal Amplitude = %0.2f\n',sig_amp);
fprintf('Standard Deviation of Noise = %0.2f\n',sig);

% SNR calculations
% 1. Fix Sigma = 1 and Vary Amplitude of Signal
figure(3); clf; 
A = linspace(1,100,15); %Amplitudes
S = 0:1:8; %Standard Deviations
ncolor = length(S-1); % |S| - 1 colors
colors = hsv(ncolor);
styles = {'+-','o-','*-','.-','x-','s-','d-','^-','v-','>-','<-','p-','h-'};
STDString = {};
for k = 2:length(S) %start from 2 b/c don't need sigma = 0
    sig1_mat = zeros(1,num) + S(k);
    sig2_mat = sig1_mat;
    u1_mat = zeros(1,num);
    u2_mat = u1_mat;
    nor1 = normrnd(u1_mat,sig1_mat);
    nor2 = normrnd(u2_mat,sig2_mat);
    ray_made = sqrt(nor1.^2 + nor2.^2);
    inphase = nor1;
    quad = nor2;
    Spower = [];
    Npower = mean(ray_made.^2);
    SNR = [];
    
    %snr
    for m = 1:length(A)
        i = 1;
        signal = zeros(1,length(inphase));
        for j = 1:length(t)
            if t(j) > Delay(i)
                signal(j) = signal(j) + A(m);
                i = i + 1;
                if i > length(Delay)
                    break;
                end
            end
        end
        signalPower = mean(signal.^2);
        noisePower = mean(inphase.^2 + quad.^2);
        SNR = [SNR signalPower/noisePower];
    end
    SNR = 10*log10(SNR);
    h = plot(A,SNR,styles{mod(k,numel(styles))+1},'Color',colors(k,:)); hold on;
    STDString{end+1} = ['\sigma',sprintf(' = %0.2f',S(k))];
end
xlabel('Amplitude'); ylabel('SNR(dB)'); grid minor;
title('Amplitude vs SNR'); legend(STDString);
