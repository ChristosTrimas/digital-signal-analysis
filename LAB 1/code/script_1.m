%lab1
%Christos Trimas 2016030054
%Kyriakos Christodoulidis 2016030025

clc;
clear all;;
close all;

%A
nx = -15:15; %time of X
X  = 4*(nx>=4) + (-nx+8>=0).*(-nx+8); %x[n]=u[n-4]+r[-n+8]
figure;
subplot(2,1,1)
stem(nx,X);
xlabel('n');
ylabel('x[n]');
title('x[n] signal');

nh = -15:10; %time of h
h = (1/2).^abs(nh); %h[n]=(1/2)^|n|
subplot(2,1,2);
stem(nh,h);
xlabel('n');
ylabel('h[n]');
title('h[n] signal');

%zero padding
ny = nx(1) + nh(1):nx(end) + nh(end); %time of conv

%lengths
LenY = length(ny); %length of convolution
LenX = length(X); %length of X
LenH = length(h); %length of h

%filling with zeros X so we can do the multiplications
X0 = [zeros(1,LenH-1) X zeros(1,LenH-1)]; %zero padding
% wide vector of time for X0
nx0 = nx(1)-(LenH-1):nx(end) + (LenH-1);

figure;
subplot(2,1,1);
stem(nx,X)
xlabel('n');
ylabel('x[n]');
title('x[n] signal]');

subplot(2,1,2);
stem(nx0,X0);
xlabel('n');
ylabel('x[n]');
title('x[n] signal with zero padding]');

%for the convolution i need h[-n]
h_rev = h(end:-1:1);

%same time limits as X0
nx0 = nx(1) - (LenH-1):nx(end) + (LenH - 1);

figure;
subplot(2,1,1);
stem(nx0,X0); %signal with zeros

%metakinhsh kata n
for i=1:LenY %length of convolution
    h_rev0 = [zeros(1,(i-1)) h_rev zeros(1,(LenY-i))];
%     figure;
%     subplot(2,1,2);
%     stem(nx0,h_rev0);     %those lines are if we want to see the signal
%     moved by n
    y1(i) = sum(X0.*h_rev0);
end


subplot(2,1,1);
stem(ny,y1);

%convolution using the conv command
y2 = conv(X,h);
subplot(2,1,2);
stem(ny,y2);

%B
f1 = 25;
f2 = 50;
n = -15:15;
y3 = cos(f1*n);
y4 = sin(2*pi*f2*n);

figure;
subplot(2,1,1);
stem(n,y3);
ylabel('y3[n]');
xlabel('n');

subplot(2,1,2);
stem(n,y4);
xlabel('n');
ylabel('y[n]');

tf = [n(1)+n(1) : n(end) + n(end)];
z = conv(y3,y4);

Fs = 2000;
DT = 1/Fs;
N = 1024;
F = [-Fs/2:Fs/N:Fs/2-Fs/N];

Z = fft(z,N);
Y3 = fft(y3,N);
Y4 = fft(y4,N);

Z0 = Y3.*Y4;

figure;
subplot(2,1,1);
stem(F,Z);
xlabel('Hz');
ylabel('Z(F)');
title('Fourier of conv y3,y4');

subplot(2,1,2);
stem(F,Z0);
xlabel('Hz');
ylabel('Z0(f)');
title('Y3 * Y4');

