%Exercise 3
%Christos Trimas 2016030054
%Kuriakos Christodoulidis 2016030025

clc;
clear all;
close all;

%Creation of lowpass Butterworth filter

fs = 10000;      %sampling frequency
Wp = 2*pi*3000;  % (passband frequency corner)
Ws = 2*pi*5000;  % (stopband frequency corner)
Rp = 3;          %ripple
Rs = 30;         %attenuation
[n,Wn] = buttord(Wp,Ws,Rp,Rs, 's'); %n:order of filter,Wn:3dB frequency

[z, p, k]=buttap(n);%Finds zeros and poles of the analog Butterworth filter

[num, den] = zp2tf(z, p, k);%Zero-pole to transfer function conversion.

[NUMT,DENT] = lp2lp(num,den,Wn);%transforms the lowpass filter  prototype
%with unity cutoff frequency of 1 rad/sec to a lowpass filter with cutoff
%frequency Wo
N = 2048;
f = linspace(0, fs/2, 2048);

f1 = freqs(NUMT,DENT,2*pi*f);%H(jw) analog

[num_z,den_z] = bilinear(NUMT,DENT,fs);%frequency wrapping
 
f2= freqz(num_z, den_z, 2048, fs);%H(e^jw) digital

figure
plot(f, (20*log10(abs(f1))), 'r--');
hold on;
plot(f, (20*log10(abs(f2))), 'b.');
legend('Analog filter', 'Digital filter');
xlabel('Frequency (Hz)');
ylabel(' |H(jù)|,  |H(e^j^ù)|');
title('Butterworth Lowpass Filter');
hold off;

Rs2 = 50;
[n2,Wn2] = buttord(Wp,Ws,Rp,Rs2, 's'); %n:order of filter,Wn:3dB frequency
Wn2/100
[z2, p2, k2]=buttap(n2);%Finds zeros,poles of the analog Butterworth filter

[num2, den2] = zp2tf(z2, p2, k2);%Zero-pole to transfer function conversion.

[NUMT2,DENT2] = lp2lp(num2,den2,Wn2);
%transforms the lowpass filter  prototype
%with unity cutoff frequency of 1 rad/sec to a lowpass filter with cutoff
%frequency Wo
N = 2048;
f = linspace(0, fs/2, 2048);

f1_2 = freqs(NUMT2,DENT2,2*pi*f);%H(jw) 

[num_z2,den_z2] = bilinear(NUMT2,DENT2,fs);%Transform into z 
 
f2_2= freqz(num_z2, den_z2, 2048, fs);%H(e^jw) 

figure
plot(f, 20*log10(abs(f1_2')), 'r--');
hold on;
plot(f, 20*log10(abs(f2_2')), 'c.');
legend('Analog filter', 'Digital filter');
xlabel('Frequency (Hz)');
ylabel(' |H(jù)|,  |H(e^j^ù)|');
title('Butterworth Lowpass Filter');
hold off;

%2nd exercise

%Creation of highpass Chebyshev filter
n1 = 2;  %filter order(1)
n2 = 16; %filter order(2)

Ts1=0.2;  %sampling period
fs1=1/Ts1; %sampling frequency
Wc=2;    %cutoff frequency
fc = 1/pi;
Rp1=3;    %passband ripple
N1=256;   %samples
Wp1=Wc/(fs1*pi); %cutoff frequency normalized 
w=0:1/(N1-1):1;    %axis of the frequency

[b_1,a_1] = cheby1(n1,Rp1,Wp1,'high');%highpass filter with order n1
f1_dig=freqz(b_1,a_1,N1);

[b_2,a_2] = cheby1(n2,Rp1,Wp1,'high');%highpass filter with order n2
f2_dig=freqz(b_2,a_2,N1);

figure
plot(w, 20*log10(abs(f1_dig)), 'r--');
hold on;
plot(w, 20*log10(abs(f2_dig)), 'b');
legend('n=2', 'n=16');
xlabel('Frequency (radians/sample)');
ylabel(' |H(jù)|(dB),  |H(e^j^ù)|(dB)');
title('Chebyshev highpass filter');
hold off;

%3rd exercise
fs=10000; %sampling frequency
Ts=1/fs;
f1=500/pi;
f2=8000/pi;
f3=15000/pi;
N2=500; %samples
n=0:N2-1;
x1=1+cos(2*pi*f1*n*Ts)++cos(2*pi*f2*n*Ts)++cos(2*pi*f3*n*Ts);

figure
subplot(2,1,1);
plot(n,x1);
xlabel('n');
ylabel('x1');
title('Samples of x(t)=1+cos(1000t)+cos(16000t)+cos(30000t)');

f_axis1=-fs/2:fs/N2:fs/2-fs/N2; % axis of the frequency
x_f=fftshift(fft(x1))*Ts;    %Fourier Transform
subplot(2,1,2);
plot(f_axis1,abs(x_f));
xlabel('Frequency (Hz)');
ylabel('x_f');
title('Fourier Transform of x(t)=1+cos(1000t)+cos(16000t)+cos(30000t)');

x_filtered=filter(num_z,den_z,x1); 
figure
subplot(2,1,1);
plot(n,x_filtered);
xlabel('n');
ylabel('x_filtered');
title('Filtered 30db x(t)=1+cos(1000t)+cos(16000t)+cos(30000t)');


x_f_filtered=fftshift(fft(x_filtered))*Ts; 
subplot(2,1,2);
plot(f_axis1,abs(x_f_filtered));
xlabel('Frequency (Hz)');
ylabel('x_filtered');
title('Filter 30db Fourier Transform x(t)=1+cos(1000t)+cos(16000t)+cos(30000t)');

%3b
Ts1=0.2;
fs1=1/Ts1; %sampling frequency
fa=0.75/pi;
fb=2.5/pi;
N2=500;
n=0:N2-1;
f_axis2=-fs1/2:fs1/N2:fs1/2-fs1/N2;
x2=1+cos(2*pi*fa*n*Ts1)+cos(2*pi*fb*n*Ts1);

figure
subplot(2,1,1);
plot(n,x2);
xlabel('n');
ylabel('x2');
title('Samples of x(t)=1+cos(1.5t)+cos(5t)');

x_fn=fftshift(fft(x2))*Ts1; %Fourier Transform
subplot(2,1,2);
plot(f_axis2,abs(x_fn));
xlabel('Frequency (Hz)');
ylabel('x_fn');
title('Fourier Transform of x(t)=1+cos(1.5t)+cos(5t)');

figure
x2n_f=filter(b_2,a_2,x2); 
subplot(2,1,1);
plot(n,x2n_f);
xlabel('n');
ylabel('x2n_f');
title('Filtered x(t)=1+cos(1.5t)+cos(5t)');

x_f2n_f=fftshift(fft(x2n_f))*Ts1;  
subplot(2,1,2);
plot(f_axis2,abs(x_f2n_f));
xlabel('Frequency (Hz)');
ylabel('x_f2_f');
title('Filtered  Fourier Transform of x(t)=1+cos(1.5t)+cos(5t)');
