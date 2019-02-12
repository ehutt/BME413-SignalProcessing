
% Question 1 

clear all 
close all
% a) 
f = 4; %4 Hz frequency  
fs = 1000; %sampling frequency 
N = 1000; %num sampling points
t = (0:N-1)/fs; %time vector 
sine = sin(2*pi*f*t); 
sq = square(2*pi*f*t); 
cosine = cos(2*pi*f*t); 
figure; 
plot(sine, 'b'); hold on;
plot(cosine, 'r'); hold on; 
plot(sq,'k')
xlabel('Time (seconds)') 
ylabel('Signal') 
ylim([-1.25 1.25])

% b)  
r_sin = mean(sine.*sq); 
r_cos = mean(cosine.*sq); 

% C)  
%manual calculation 
% r_sin = (1/T)*integral(sine*sq) over [0,T] T = 1/4 = 250ms 
%       = integral(sine*sq) over [0,.5T] + integral(sine*sq) over [0.5T,T]
%       = -cos(0) + cos(0.5*T) + cos(0.5*T) - cos(T) 
%       = -1 +1 +1 -1
%       = 0

% r_cos = (1/T)*integral(cosine*sq) over [0,T] T = 1/4 = 250ms 
%       = integral(cos*sq) over [0,.5T] + integral(cos*sq) over [0.5T,T]
%       = sin(0) - sin(0.5*T) - sin(0.5*T) + sin(T) 
%       = sin(0) - sin(T)
%       = 0 

% Question 2 

load dataforEEG.mat 
N = length(eeg); %num samples
interval = 16; %16 sec signal 

% a)
%sampling frequency: 
fs = N/interval; 
t = (1:length(eeg))/fs;

% b)

rmax = (1:20); 
for i=1:20
    x = sin(2*pi*i*t); 
    r = xcorr(eeg,x); 
    rmax(i) = max(r); 
end
figure;
plot(rmax)
xlabel('Frequency (Hz)') 
ylabel('Max Xcorr') 

%max peak is in 7Hz range, this is the alpha wave 

% Question 3  
load dataforResp 

% a) 
N = length(resp); %num samples
fs = N / 125; %sampling frequency for 125s interval 
f = (1:N)*fs/N;%frequency vector 
X = fft(resp); %complex fourier transform 

%only plot spectrum up to 2Hz 
index = round(2*(N/fs));  

%plot magnitude 
subplot(2,1,1);
plot(f(1:index-1),abs(X(2:index)),'k'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');

%unwrap smooths out discontinuities 
phase = unwrap(angle(X)); 

%plot phase 
subplot(2,1,2); 
plot(f(1:index-1),phase(2:index),'k'); 
xlabel('Frequency (Hz)'); 
ylabel('Phase'); 

% b) 

%find maximal peak and its frequency/time location
[peak, n_peak] = max(abs(X(2:index))); 
max_freq = f(n_peak); %breathing frequency 
max_time = 1/max_freq; 
breaths_per_minute = 60 / max_time; 

