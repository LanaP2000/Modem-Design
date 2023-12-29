clear all;
close all;
clc;

fc = 100; % Carrier frequency at the transmitter side
symbTime = 0.01; % Symbol time = 0.01s
M = 10;
Ts = 1 / (M * fc);

sampPerSymb = symbTime / Ts;
t = 0:Ts:symbTime-Ts;

EbN0_dB = 2:1:8;
EbN0 = 10.^(EbN0_dB/10);
Eb = 1;
BER = zeros(1,length(EbN0_dB));
beta = 0.5;
filter_length = 10;

ps = srrc(filter_length, beta, sampPerSymb);  % Squared root raised cosine filter

lp = fir1(48, 0.1);

filter_delay = sampPerSymb * 10 * 2 + 1;

L = 10000;
N_symbols = L / 2; % For QPSK

for i = 1:length(EbN0_dB)
    EbN0_dB(i)
    N0 = 1 / EbN0(i); % With Eb = 1, N0 = 1/(EbN0)
    error = 0;
    noOfL = 0;
    a = sqrt(Eb);
    while (error < 500) && (noOfL < 1e4)
        noOfL = noOfL + 1;
        bits = randi([0, 1], 1, L); % Random bit sequence of length L

        bits_I = bits(1:2:end);
        bits_Q = bits(2:2:end);

        symbols_I = a * (bits_I * 2 - 1);
        symbols_Q = a * (bits_Q * 2 - 1);

        symbols_I_up = upsample(symbols_I, sampPerSymb);
        symbols_Q_up = upsample(symbols_Q, sampPerSymb);

        signal_I_shaped = conv(ps, symbols_I_up); % Convolve SRRC shape with data
        signal_Q_shaped = conv(ps, symbols_Q_up);

        signal_I = signal_I_shaped .* cos(2 * pi * fc * [0:length(signal_I_shaped)-1] * Ts) * sqrt(2);
        signal_Q = signal_Q_shaped .* sin(2 * pi * fc * [0:length(signal_I_shaped)-1] * Ts) * sqrt(2);

        tx = signal_I + signal_Q; % Transmitted signal

        % Observed signal at the receiver, oversampled by a factor of 100
        [rx_I, rx_Q, freq_offset, phase_offset] = rx_offset(tx, fc, Ts, N0, 1, 0);

        % SRRC at the receiver
        rx_I_shaped = conv(ps, rx_I);
        rx_Q_shaped = conv(ps, rx_Q);

        rp = [rx_I_shaped, rx_Q_shaped];

        mu1 = .01; 
        mu2 = .003; % algorithm step sizes
        f0 = 101; % assumed freq at receiver
        lent = length(t); 
        th1 = zeros(1, lent); % initialize estimates
        th2 = zeros(1, lent); 
        carest = zeros(1, lent);
        for k = 1 : lent - 1 % combine top PLL th1
            th1(k + 1) = th1(k) - mu1 * rp(k) * sin(4 * pi * f0 * t(k) + 2 * th1(k));
        end

        remd = exp(-1i * th1);

        % Remove frequency offset
        rx_I_shaped_ = rx_I_shaped .* remd(end);
        rx_Q_shaped_ = rx_Q_shaped .* remd(end);

        % Downsample to get to the sample
        rx_symbols_I = rx_I_shaped_(filter_delay:sampPerSymb:filter_delay+sampPerSymb*N_symbols-1);
        rx_symbols_Q = rx_Q_shaped_(filter_delay:sampPerSymb:filter_delay+sampPerSymb*N_symbols-1);

        bits_I_demod = double(rx_symbols_I > 0);
        bits_Q_demod = double(rx_symbols_Q > 0);

        bits_demod = [bits_I_demod; bits_Q_demod];
        bits_demod = reshape(bits_demod, 1, numel(bits_demod));

        error = error + sum(bits(L/5+1:end) ~= bits_demod(L/5+1:end));
    end
    BER(i) = error / L / noOfL / 0.8;
end

BER
figure;
scatter(real(rx_symbols_I), real(rx_symbols_Q))
figure;
semilogy(EbN0_dB, BER, 'LineWidth', 1.5)
BER


