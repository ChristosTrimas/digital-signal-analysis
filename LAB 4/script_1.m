%Lab 4
%Exercise 2
%Christos Trimas 2016030054
%Kuriakos Christodoulidis 2016030025

clc;
clear all;
close all;

%1
f_axis = @(w,Fs) 0:Fs/(2*length(w)):Fs/2-Fs/(2*length(w));

Wc = 0.4*pi;   
%cutoff frequency
Fc = Wc/(2*pi);
%sampling frequency
Fs = 100;
%window
N = 21;
%normalized cutoff frequency
Wn = Fc/(Fs/2);

%creation of the filters
hammFilter = fir1(N-1,Wn,hamming(N));
rectFilter = fir1(N-1,Wn,rectwin(N));

[h1,w1] = freqz(hammFilter,N);
[h2,w2] = freqz(rectFilter,N);      

hamming_freq = f_axis(w1,Fs);
rect_freq = f_axis(w2,Fs);

figure;
plot(hamming_freq,abs(h1),rect_freq,abs(h2));
xlabel('F(Hz)');
ylabel('Magnitude');
legend('Hamming','Rectangular');
title('Frequency Response of both Filters');

