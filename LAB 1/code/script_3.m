clc;
clear all;
close all;

DT = 0.01;
mint = 0;
maxt = 0.5;
t = mint:DT:maxt;
Fs = 200;
Ts = 1/Fs;

x = 10*cos(2*pi*20*t) - 4*sin(2*pi*40*t);

figure;
subplot(3,1,1);
plot(t,x);
xlabel('t');
ylabel('x(t)');
title('Signal x(t)');

minN = ceil(mint/Ts);
maxN = floor(maxt/Ts);
n = minN:maxN;

Xs = 10*cos(2*pi*20*n*Ts) - 4*sin(2*pi*40*n*Ts);

subplot(3,1,2);
stem(n*Ts,Xs);
title('Signal Xs(N=128)');
xlabel('n');
ylabel('Xs');

N = 128;
X = fftshift(fft(Xs,N)*Ts);
F = -Fs/2:Fs/N:Fs/2-Fs/N;

subplot(3,1,3);
plot(F,abs(X));
title('Sima sth sixnothta');
xlabel('Hz');
ylabel('|X(F)|');

% %B
Fs1 = 8000;
Ts1 = 1/Fs1;
N1 = 2048;

F1 = -Fs1/2:Fs1/N1:Fs1/2-Fs1/N1;
i = 1;
figure;

for f=100:125:475
    x = sin(2*pi*f/Fs1*n + pi/2);
    X = fftshift(fft(x,N1)*Ts1);
    
    subplot(4,1,i);
    i = i +1;
    plot(F1,abs(X));
    title(['Fasmata Gia f = ' num2str(f)]);
    xlabel('Hz');
    ylabel('|X(F)|');
end;

i = 1;
figure;

for f=7525:125:7900
    x = sin(2*pi*f/Fs1*n + pi/2);
    X = fftshift(fft(x,N1)*Ts1);
    
    subplot(4,1,i);
    i = i + 1;
    plot(F1,abs(X));
    title(['Fasmata Gia f = ' num2str(f)]);
    xlabel('Hz');
    ylabel('|X(F)|');
end;
