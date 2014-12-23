
%% Problem 3.  BPSK Signal
format compact
clc
clear all
% Request Input of n, fs, f0, fb, and Eb/N0
fs=input(' Sample frequency (Hz)?, R: ');  %Sampling Freq
f0=input(' Digital carrier frequency(Hz)?, |f0| <= 1/2, R: '); %Carrier Freq
fb=input(' Digital data frequency(Hz)?, |fb| <= |f0| <= 1/2, R: '); %Data Freq
EbN0=input(' Eb/N0?, ratio, R: ');
n=input(' Number of samples?, n > 1, I: ');
sd=1; % Standard deviation of 1
A=2*sd*sqrt(fb*EbN0);
t=0:n-1;  %Time variable in seconds


%%	Form the carrier

x = A*cos(2*pi*f0*t);

% Produce a random sequence of {0,1}
% rand('uniform');
j = 1;
for i=1:n;
 if fix((i-1)*fb) < fix(i*fb);
  r=rand ;
  if r <0.5;
     j=-1;
   else; 
     j=1;
  end;
 end;
 p(i)=j;
end;

y = p.*x;


Eb=(10*log10((A^2/2)/(fb*fs)))*ones(size(y));
N0=(10*log10(2*sd^2/fs))*ones(size(y)); % SD of randn is 1
figure(1);
hold on;

axis([0 500 -70 0]);
ZZ=(1/fs)*fft(z,n);
ZD=10*log10(2*(fs/n)*(abs(ZZ).^2));

f=0:fs/n:(n/2-1)*fs/n;
plot(f,ZD(1:n/2),'y');
plot(Eb,'r');
plot(N0,'g');

grid on;

 h = spectrum.welch;    % Create a Welch spectral estimator. 
 Hpsd = psd(h,z,'fs',fs);             % Calculate the PSD 
             plot(Hpsd)  % Plot the PSD.
title('Single-sided Periodogram and Welch PSD in dB/Hz versus frequency' );
hold off;


% Generate n samples of of sampled BPSK
% Carrier frequency = f0*fs
% Data Frequency = fb*fs and radom w/ equally likely 1 and 0

% Add White Gaussian Noise 


%%  test



