clc;
clear all;
close all;

DT = 0.01;
t = (0:DT:0.5) %500ms = 0,5 s

x = 5*cos(24*pi*t) - 2*sin(1.5*pi*t);

maxT = 0.5;
minT = 0;

T = [1/48 1/24 1/12];

figure;

for i = 1:3
    maxN = floor(maxT/T(i));
    minN = ceil(minT/T(i));
    
    n = minN:maxN;
    
    Xs = 5*cos(24*pi*n*T(i)) - 2*sin(1.5*pi*n*T(i));
    
    subplot(3,1,i);
    hold on;
    plot(t,x,'y');
    stem(n*T(i),Xs,'r');
    xlabel('t,b');
    ylabel('x(t),Xs[n]');
    
end;