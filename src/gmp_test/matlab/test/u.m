clc;
close all;
clear;


x = 0:0.001:0.1;
y = 1 - 1 ./ (1 + exp(100*(x-0.05)) );
figure;
plot(x,y, 'LineWidth',2);
return;