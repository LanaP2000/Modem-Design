clear all;
close all;
clc;

Ts = 1/10000; 
time = 1; 
t = Ts :Ts : time; % time vector
fc = 100; 
phoff = -0.8; % carrier freq and phase
rp = cos(4 * pi * fc * t + 2 * phoff) * sqrt(2); % simplified rec'd signal
fl = 100; 
ff = [0 .01 .05 1]; 
fa = [1 1 0 0];
h = firpm(fl, ff, fa); % LPF design
mu = .001; % algorithm step size
f0 = 100; % freq at receiver
theta = zeros(1, length(t)); 
theta(1) = 0; % initialize estimates
z = zeros(1, fl + 1); % initialize LPF
for k = 1 : length(t) - 1 % z contains past inputs
    z = [z(2 : fl + 1), rp(k) * sin(4 * pi * f0 * t(k) + 2 * theta(k))];
    update = fliplr(h) * z'; % new output of LPF
    theta(k + 1) = theta(k) - mu*update; % algorithm update
end

plot(t, theta)
