function [rx_I, rx_Q, freq_offset, phase_offset] = rx_offset(tx, fc, Ts, N0, freq_opt, phase_opt)
freq_offset = 0;
phase_offset = 0;
if freq_opt == 1
    freq_offset = randi([1,10],1,1)*1e-6*fc;
end

if phase_opt == 1    
    phase_offset = rand*2*pi - pi;
end

fc_rx = fc + freq_offset;
noise = sqrt(N0/2) * randn(1, length(tx));
rx = tx + noise;
rx_I = rx.*cos(2*pi*fc_rx*[0:length(rx)-1]*Ts + phase_offset)*sqrt(2);
rx_Q = rx.*sin(2*pi*fc_rx*[0:length(rx)-1]*Ts + phase_offset)*sqrt(2);