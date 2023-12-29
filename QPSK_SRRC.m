clearvars;
fc = 1e2; % carrier frequency at the transmitter side
symbTime = 0.1; % symbol time = 0.1s
M = 10;
Ts = 1/M/fc; %

sampPerSymb = symbTime/Ts;
t = 0:Ts:symbTime-Ts;

EbN0_dB = 2:1:8;
EbN0 = 10.^(EbN0_dB/10);
Eb = 1;
BER = zeros(1,length(EbN0_dB));
beta = 0.5;
filter_length = 10;

ps=srrc(filter_length,beta,sampPerSymb);  % squared root raised cosine filter

lp = fir1(48,0.1);

filter_delay = sampPerSymb*10*2 + 1;

L = 10000;
N_symbols = L/2; % for QPSK
for i=1:length(EbN0_dB)
    EbN0_dB(i)
    N0 = 1/EbN0(i);              % with Eb = 1, N0 = 1/(EbN0)
    error = 0;
    noOfL = 0;
    a = sqrt(Eb);
    while (error < 500) && (noOfL<1e4)
        noOfL = noOfL + 1;
        bits = randi([0,1],1,L); % random bit sequence of length L

        bits_I = bits(1:2:end);
        bits_Q = bits(2:2:end);
        
        symbols_I = a*(bits_I*2 - 1);
        symbols_Q = a*(bits_Q*2 - 1);
        
        symbols_I_up = upsample(symbols_I,sampPerSymb);
        symbols_Q_up = upsample(symbols_Q,sampPerSymb);
                        
        signal_I_shaped = conv(ps,symbols_I_up);     % convolve SRRC shape with data
        signal_Q_shaped = conv(ps,symbols_Q_up);
         
        signal_I = signal_I_shaped.*cos(2*pi*fc*[0:length(signal_I_shaped)-1]*Ts)*sqrt(2);
        signal_Q = signal_Q_shaped.*sin(2*pi*fc*[0:length(signal_I_shaped)-1]*Ts)*sqrt(2);
        
        tx = signal_I + signal_Q; % transmitted signal
        
        % Your tasks are here
        % rx_offset(tx, fc, Ts, N0, freq_offset_opt, phase_offset_opt);
        % set freq_offset_opt to 1: to introduce a random frequency offet
        % set phase_offset_opt to 1: to introduce a random phase offset.

        % observed signal at the received, oversampled by a factor of 100
        [rx_I, rx_Q] = rx_offset(tx, fc, Ts, N0, 0, 0); 

        % SRRC at the receiver
        rx_I_shaped = conv(ps,rx_I);
        rx_Q_shaped = conv(ps,rx_Q);
        
        % downsampling to get to the sample
        rx_symbols_I = rx_I_shaped(filter_delay:sampPerSymb:filter_delay+sampPerSymb*N_symbols-1);
        rx_symbols_Q = rx_Q_shaped(filter_delay:sampPerSymb:filter_delay+sampPerSymb*N_symbols-1);
       
        % uncomment to scatter plot to view the constellation if there is
        % freqeueny offset and/or phase offset.
%         scatter(rx_symbols_I, rx_symbols_Q);

        bits_I_demod = double(rx_symbols_I > 0); %PAM4_demod_GRAY(rx_symbols_I/a);
        bits_Q_demod = double(rx_symbols_Q > 0); %PAM4_demod_GRAY(rx_symbols_Q/a);
        
        bits_demod = [bits_I_demod;bits_Q_demod];
        bits_demod = reshape(bits_demod, 1, numel(bits_demod));

        error = error + sum(bits(L/5+1:end) ~= bits_demod(L/5+1:end));
    end
    BER(i) = error/L/noOfL/0.8;
end
BER
semilogy(EbN0_dB,BER,'LineWidth',1.5)
figure;
scatter(rx_symbols_I, rx_symbols_Q);


