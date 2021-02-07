clc; clear all;
% Plot of Msg Signal
fs = 3000;
fm = 50;
ts = 1/fs;
wm = 2*pi*fm;
t = 0:ts:2;
x_signal = cos(wm*t) + 2*cos(2*wm*t) + 4*cos(3*wm*t) + 5*cos(5*wm*t); % Message Signal
figure(1);
plot(t,x_signal);
xlabel('Time (s)')
ylabel('x(t)')
title('Plot of Message Signal')
%%
% Conventional AM signal u=0.8
Ac = 1;
u = 0.8; % Modulation index
fc = 20*fm;
wc = 2*pi*fc;
x1_signal = Ac*cos(wc*t); % carrier 
% figure(2);
% plot(t,x1_signal);
x2_signal = (1+(u*x_signal)).*x1_signal;
figure(2);
plot(t,x2_signal);
xlabel('Time (s)')
ylabel('x_{AM}(t)')
title('Plot of Modulated Signal u=0.8')

% AM signal u=1.2
u = 1.2;
x3_signal = (1+(u*x_signal)).*x1_signal;
figure(3);
plot(t,x3_signal);
xlabel('Time (s)')
ylabel('x_{AM}(t)')
title('Plot of Modulated Signal u=1.2')
%%
% Spectrum of x_signal / x(t)
fft_x = fftshift(fft(x_signal));
fft_x_mag = abs(fft_x);
faxis=linspace(-fs/2,fs/2-fs/length(fft_x),length(fft_x));
figure(4);
plot(faxis,fft_x_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of x(t)')

% Spectrum of x2_signal / x_AM(t)
fft_x2 = fftshift(fft(x2_signal));
fft_x2_mag = abs(fft_x2);
fft_x2_phase = angle(fft_x2);
faxis2=linspace(-fs/2,fs/2-fs/length(fft_x2),length(fft_x2));
figure(5);
plot(faxis2,fft_x2_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of x_{AM}(t)')
%%
% DSB-SC AM Generation
% Time domain
x_dsb_signal = x_signal.*x1_signal; % x(t)*Ac*cos(wct)
figure(6);
plot(t,x_dsb_signal);
xlabel('Time (s)')
ylabel('x_{dsb}(t)')
title('Plot of x_{DSB}(t)')

% Freq domain
fft_xdsb = fftshift(fft(x_dsb_signal));
fft_xdsb_mag = abs(fft_xdsb);
faxisdsb=linspace(-fs/2,fs/2-fs/length(fft_xdsb),length(fft_xdsb));
figure(7);
plot(faxisdsb,fft_xdsb_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of x_{DSB}(t)')
%%
% Demodulation
fs = 3000;
ts = 1/fs;
t = 0:ts:2;
fc = 1000;
wc = 2*pi*fc;
x1_signal = cos(wc*t);
x_sig = cos(2*pi*50*t); % Message signal
x_dsb_sig = x_sig.*x1_signal; % Dsb modulated signal

fft_x_sig = fftshift(fft(x_sig));
fft_x_sig_mag = abs(fft_x_sig);
faxis=linspace(-fs/2,fs/2-fs/length(fft_x_sig),length(fft_x_sig));
figure(8);
plot(faxis,fft_x_sig_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Message Signal')
xlim([-60,60]);

figure(9);
plot(t,abs(x_sig));
xlabel('Time (s)')
ylabel('Magnitude')
title('Plot of Message Signal')

% Csignal cos(2*pi*50t)
c_sig1 = cos(2*pi*50*t); % carrier signal
x_c_sig1 = x_dsb_sig.*c_sig1;
fft_sig1 = fftshift(fft(x_c_sig1));
[B,A,k] = butter(1,0.03,"low");
fft_y1 = filter(B,A,fft_sig1); % Pass through low-pass filter
fft_y1_mag = abs(fft_y1);
figure(10);
plot(faxis,fft_y1_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Demodulated Signal for carrier cos(2π50t)')

y1 = ifft(fft_y1);
figure(11);
plot(t,abs(y1));
xlabel('Time (s)')
ylabel('Magnitude')
title('Plot of Demodulated Signal for carrier cos(2π50t)')

% Csignal cos(2*pi*50t + pi/4)
c_sig2 = cos(2*pi*50*t + pi/4);
x_c_sig2 = x_dsb_sig.*c_sig2;
fft_sig2 = fftshift(fft(x_c_sig2));
[B,A,k] = butter(1,0.03,"low");
fft_y2 = filter(B,A,fft_sig2);
fft_y2_mag = abs(fft_y2);
figure(12);
plot(faxis,fft_y2_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Demodulated Signal for carrier cos(2π50t+π/4)')
ylim([0,1400]);

y2 = ifft(fft_y2);
figure(13);
plot(t,abs(y2));
xlabel('Time (s)')
ylabel('Magnitude')
title('Plot of Demodulated Signal for carrier cos(2π50t+π/4)')
ylim([0,1.2]);

% Csignal cos(2*pi*50t + pi/2)
c_sig3 = cos(2*pi*50*t + pi/2);
x_c_sig3 = x_dsb_sig.*c_sig3;
fft_sig3 = fftshift(fft(x_c_sig3));
[B,A,k] = butter(1,0.03,"low");
fft_y3 = filter(B,A,fft_sig3);
fft_y3_mag = abs(fft_y3);
figure(14);
plot(faxis,fft_y3_mag);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of Demodulated Signal for carrier cos(2π50t+π/2)')
ylim([0,1400]);

y3 = ifft(fft_y3);
figure(15);
plot(t,abs(y3));
xlabel('Time (s)')
ylabel('Magnitude')
title('Plot of Demodulated Signal for carrier cos(2π50t+π/2)')
ylim([0,1.2]);