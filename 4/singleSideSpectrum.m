function [P,f] = singleSideSpectrum(signal)
Fs = 10;
L = length(signal);
Y = fft(signal);
P2 = abs(Y/L);
P = P2(1:L/2+1);
P(2:end-1) = 2 * P(2:end-1);
f = Fs * (0:(L/2))/L;
end