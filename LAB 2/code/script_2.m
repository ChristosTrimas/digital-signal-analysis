%Lab 2 
%Exercise 2
%Christos Trimas 2016030054
%Kuriakos Christodoulidis 2016030025

clc;
clear all;
close all;

X = [4 -3.5 0];
Y = [1 -2.5 1];

[A,B,C] = residuez(X,Y); %evresh polwn/sintelestwn

syms z; %orismos tou migadikou z

H_1 = A(1)/(1-B(1)*z^(-1));
H_2 = A(2)/(1-B(2)*z^(-1));
H = H_1 + H_2;
pretty(H)

Hn = iztrans(H)




