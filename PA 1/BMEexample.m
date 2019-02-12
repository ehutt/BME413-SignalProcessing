%BME discussion example 


%%%%%% Q1 
clear all 
close all 
load dataforResp 
N = length(resp); %num samples
fs = 100; %sample frequency = 100hz 
f = (1:N)*fs/N;%frequency vector 

plot(t,resp,'k'); %plot signal in time
X = fft(resp); %complex fourier transform 
m_plot = round(4*(fs/N)); %only show below 4hz frequency


%plot first 4hz 
%start at X(2) bc X(1) is DC component 
%peaks of this plot are the first few harmonics 
subplot(2,1,1); 
plot(f(1:m_plot-1),abs(X(2:m_plot)),'k'); 
xlabel('Frequency (Hz)'); 

%unwrap smooths out discontinuities 
phase = unwrap(angle(X)); 

%plot phase 
subplot(2,1,2); 
plot(f(1:m_plot-1),phase(2:m_plot),'k'); 
xlabel('Frequency (Hz)'); 

%only care about first harmonic, tallest peak 
%find maximal peak and its frequency/time location
[peak, m_peak] = max(abs(X(2:m_plot))); 
max_freq = f(m_peak); %breathing frequency 
max_time = 1/max_freq; 


%%%%%%% Q2 

%alpha wave 

clear all 
close all 

load dataforEEG 
fs = 30; %sampling frequency 30 mega hz 
t = (1:length(eeg))/fs; 

for i=1:15 
    f(i) = i; 
    x = cos(2*pif(i)*t); 
    r = xcorr(eeg,x); 
    rmax(i) = max(r); 
end



