%Lab 2 
%Exercise 2
%Christos Trimas 2016030054
%Kuriakos Christodoulidis 2016030025

clc;
clear all;
close all;

%b
fs=1;
Ts=1/fs;
num1=[0.2 0];
den1=[1 -0.9];
G1=tf(num1,den1,Ts);

num2=[0 1];
den2=[1 0.2];
G2=tf(num2,den2,Ts);


H1=G1*G2
zero = 0;
poles = [0.9 ; -0.2 ];
figure
zplane(zero,poles)
title('Pole/Zero Plot for H(z)=0.2z/z^{2}-0.7z-0.18');

%d
a=[-pi:pi/128:pi]; %diasthma 
numerator=[0 0.2 0];
denominator=[-0.18 -0.7 1];
figure
freqz(numerator,denominator,a)
figure
freqz(numerator,denominator)

%e
num3=[0 1];
den3=[1 -1];
G3=tf(num3,den3,Ts);
H2=H1*G3
denominator2=[0.18 0.52 -1.7 1 ];
figure
freqz(numerator,denominator2,a)


%2
X = [4 -3.5 0];
Y = [1 -2.5 1];

[A,B,C] = residuez(X,Y); %evresh polwn/sintelestwn

syms z; %orismos tou migadikou z

H_1 = A(1)/(1-B(1)*z^(-1));
H_2 = A(2)/(1-B(2)*z^(-1));
H = H_1 + H_2;
pretty(H)

Hn = iztrans(H)
