function [spec,f]=plotfft(trace,dt)
NL=length(trace);
NFFT=2^nextpow2(NL);
spec=fft(trace,NFFT)/NL;
f=(1/(2*dt))*(linspace(0,1,NFFT/2+1));
spec=spec(1:NFFT/2+1,:);