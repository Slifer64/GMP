clc;
close all;
clear;

Q = 2*rand(4,1) - 1;
Q = Q / norm(Q);

omega = 2*rand(3,1) - 1;

omega_dot = 2*rand(3,1) - 1;


qDot = gmp_.rotVel_to_qLogDot(omega, Q);

omega2 = gmp_.qLogDot_to_rotVel(qDot, Q);

norm(omega - omega2)



qddot = gmp_.rotAccel_to_qLogDDot(omega_dot, omega, Q);

omega_dot2 = gmp_.qLogDDot_to_rotAccel(qddot, omega, Q);

norm(omega_dot2 - omega_dot)

% norm(omega) / norm(qDot)
