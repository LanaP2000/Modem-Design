clear all;
close all;
clc;

% dualplls.m estimation of carrier via dual-loop structure

Ts = 1/10000; 
time = 5; 
t = 0 : Ts : time - Ts; % time vector
fc = 1000; 
phoff = -2; % carrier freq and phase
rp = cos(4 * pi * fc * t + 2 * phoff); % preprocessed signal rBPF
mu1 = .01; 
mu2 = .003; % algorithm step sizes
f0 = 1001; % assumed freq at receiver
lent = length(t); 
th1 = zeros(1, lent); % initialize estimates
th2 = zeros(1, lent); 
carest = zeros(1, lent);
for k = 1 : lent - 1 % combine top PLL th1
    th1(k + 1) = th1(k) - mu1 * rp(k) * sin(4 * pi * f0 * t(k) + 2 * th1(k));
    th2(k + 1) = th2(k) - mu2 * rp(k) * sin(4 * pi * f0 * t(k) + 2 * th1(k) + 2 * th2(k));
    % with bottom PLL th2 to form estimate of preprocessed signal
    carest(k) = cos(4 * pi* f0 * t(k) + 2 * th1(k) + 2 * th2(k));
end

figure;
plot(t, th1)
figure;
plot(t, th2)
figure;
plot(t, carest)


