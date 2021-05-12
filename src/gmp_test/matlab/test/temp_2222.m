clc;
close all;
clear;

rng(0);

I3 = eye(3,3);

Jo = 5*I3 + 0.1*rand(3,3);
mo = 10;

r = rand(3,1);

Ji = Jo + mo*(norm(r)^2*I3 - r*r');

cross_r = math_.vector2ssMatrix(r);

Ji2 = Jo - mo*cross_r*cross_r;

Ji
Ji2