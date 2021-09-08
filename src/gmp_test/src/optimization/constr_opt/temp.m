clc;
close all;
clear;

kron(speye(2), [1 2; 3 4])

return

n = 3;
m = 2;

Hi = blkdiag(0.2*eye(m,m), eye(n,n));
Hcell = repmat({Hi}, 1, 3);
H = blkdiag(Hcell{:});

H
return


N = 4;

n = 2; % state dim
m = 1; % control input dim

In = eye(n,n);

Ai = [1 2; 3 4];
Bi = [0; 0.5];

Aeq = zeros( n*N, (n+m)*N );
beq = zeros( n*N , 1);

Aeq(1:n, 1:m) = -Bi;
Aeq(1:n, m+1:m+n) = In;
k = m+1;
for i=1:N-1
    i1 = 1 + i*n;
    i2 = i1 + (n-1);

    Aeq(i1:i2, k:k+n-1) = -Ai;
    Aeq(i1:i2, k+n:k+n+m-1) = -Bi;
    Aeq(i1:i2, k+n+m: k+2*n+m-1) = In;

    k = k+n+m;
end

z0 = [p; p_dot];
ud_0 = K*p_ref + D*dp_ref + ddp_ref;

Beq(1:n) = Ai*z0 + Bi*ud_0;

xi = x;
xi_dot = x_dot;
xi_ddot = x_ddot;
for i=1:N-1
    i1 = 1 + i*n;
    i2 = i1 + (n-1);

    xi = xi + xi_dot*dt;
    xi_dot = xi_dot + xi_ddot*dt;

    pd_i = gmp2.getYd(xi);
    dpd_i = gmp2.getYdDot(xi, xi_dot);
    ddpd_i = gmp2.getYdDDot(xi, xi_dot, xi_ddot);
    ud_i = K*pd_i + D*dpd_i + ddpd_i;

    Beq(i1:i2) = Bi*ud_i;
end

