% plotspec(x,t,Ts) plots the spectrum of the signal x
% Ts = time (in seconds) between adjacent samples in x
% t is the time indexing fo the signal
function plotspec(x,Ts,t)
N=length(x);                               % length of the signal x
if nargin==2
    t=Ts*(1:N);                            % define a time vector
end
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);       % frequency vector
fx=fft(x(1:N));                            % do DFT/FFT
fxs=fftshift(fx);                          % shift it for plotting
subplot(2,1,1), plot(t,x)                  % plot the waveform
xlabel('seconds'); ylabel('amplitude')     % label the axes
subplot(2,1,2), plot(ssf,abs(fxs))         % plot magnitude spectrum
xlabel('frequency'); ylabel('magnitude')   % label the axes
