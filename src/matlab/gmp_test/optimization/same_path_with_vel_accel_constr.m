clc;
close all;
clear;

rng(0);

set_matlab_utils_path('../');

%% Generate data
Ts = 0.002;
[Timed, Pd_data, dPd_data, ddPd_data] = createData(0, 10, Ts, [0; 0; 0], [0.8; 0.75; 0.9], [0.25; -0.6; 0.25]);
% load('data/train_data.mat', 'Data');
% Timed = Data.Time;
% Pd_data = Data.Pos;
% dPd_data = Data.Vel;
% ddPd_data = zeros(size(dPd_data));
% Ts = Timed(2) - Timed(1);
% for i=1:size(ddPd_data,1), ddPd_data(i,:) = [diff(dPd_data(i,:)), 0] / Ts; end
p0d = Pd_data(:,1);
pgd = Pd_data(:,end);
Tf = Timed(end);

%% train GMP
N_kernels = 25;
n_dofs = size(Pd_data,1);
gmp = GMP(n_dofs, N_kernels, 2);
train_err = gmp.train('LS', Timed/Tf, Pd_data);
train_err

Time = [];
Pref_data = [];
dPref_data = [];
ddPref_data = [];
P_data = [];
dP_data = [];
ddP_data = [];
x_data = [];
dx_data = [];
ddx_data = [];

kt = 1;
ks = 2;
p0 = p0d;
pg = ks.*(pgd - p0d) + p0;

gmp.setY0(p0);
gmp.setGoal(pg);

tau = Tf / kt;

v_lim = [-0.4 0.4; -0.2 0.4; -0.4 0.4];
a_lim = [-0.15 0.15; -0.15 0.15; -0.15 0.15];

t = 0;
dt = Ts;
x = 0;
xd_dot = 1/tau;
x_dot = xd_dot;
x_dot_prev = x_dot;
x_ddot = 0;
p = p0;
dp = zeros(n_dofs,1);
ddp = zeros(n_dofs,1);

while (true)
    
%     v_ul = min( v_lim(:,2) ./ gmp.getYdDot(x,1) );
%     v_ll = gmp.getYdDot(x, 1) ./ v_lim(:,1);
%     
%     x_dot
%     v_ul
%     1e-5/(x_dot - v_ul)
%     '---------------'
%     pause

    av = 5e-3;
    fv = av*sum(1./(dp - v_lim(:,2)).^3)/3 - av*sum(1./(dp - v_lim(:,1)).^3)/3;
    fv_dot = 0;%-3*av*( 1/(y - yL_up)^4 + 1/(y - yL_low)^4 )*y_dot;
    fa = 0*1000*( sum(1./(ddp - a_lim(:,2)).^3)/3 - sum(1./(ddp - a_lim(:,1)).^3)/3 );
    
    % x_ddot = 20*(xd_dot - x_dot) + fv;
    
    x_3dot = 5*(fv - x_ddot) + 20*(xd_dot - x_dot) + fv_dot + fa;
    
    % saturate phase variable
    if (x >= 1)
       x = 1;
       x_dot = 0;
       x_ddot = 0;
    end
    
%     bi = gmp.getYdDot(x, 1);
%     H = 1;
%     f = -xd_dot;
%     A = [-bi; bi];
%     b = [-v_lim(:,1); v_lim(:,2)];
%     opt = optimoptions('quadprog','Display','off');
%     x_dot = quadprog(H,f,A,b,[],[],[],[],xd_dot, opt);
%     x_dot

    % get reference trajectory
    p_ref = gmp.getYd(x);
    dp_ref = gmp.getYdDot(x, x_dot);
    ddp_ref = gmp.getYdDDot(x, x_dot, x_ddot);
    
    % logging
    Time = [Time t];
    Pref_data = [Pref_data p_ref];
    dPref_data = [dPref_data dp_ref];
    ddPref_data = [ddPref_data ddp_ref];
    P_data = [P_data p];
    dP_data = [dP_data dp];
    ddP_data = [ddP_data ddp];
    x_data = [x_data x];
    dx_data = [dx_data x_dot];
    ddx_data = [ddx_data x_ddot];
    
    % stopping criteria
    if (norm(p-pg)<1e-3 && x>=1), break; end

    % dynamical equations
    ddp = ddp_ref + 30*(dp_ref - dp) + 200*(p_ref - p);
    
    % integration
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    x_ddot = x_ddot + x_3dot*dt;
    p = p + dp*dt;
    dp = dp + ddp*dt;
    
end

[P2_data, dP2_data, ddP2_data] = getOutput(gmp, x_data, xd_dot, 0);


%% Plot results

figure; hold on;
plot3(Pref_data(1,:), Pref_data(2,:), Pref_data(3,:), 'LineWidth',2, 'Color','green');
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
plot3(P2_data(1,:), P2_data(2,:), P2_data(3,:), 'LineWidth',2, 'Color','blue', 'LineStyle','--');
scatter3(p0(1),p0(2),p0(3), 'LineWidth',3, 'SizeData',100, 'Marker','o', 'MarkerEdgeColor','cyan')
scatter3(pg(1),pg(2),pg(3), 'LineWidth',3, 'SizeData',100, 'Marker','x', 'MarkerEdgeColor','red')
legend({'$p_{ref}$','$p$','$p_2$'}, 'interpreter','latex', 'fontsize',15, 'Position',[0.7918 0.7283 0.1494 0.1833]);
title('3D path', 'interpreter','latex', 'fontsize',20);
xlabel('$x$', 'interpreter','latex', 'fontsize',16);
ylabel('$y$', 'interpreter','latex', 'fontsize',16);
zlabel('$z$', 'interpreter','latex', 'fontsize',16);
view([-15 25.4]);
axis tight;

fig = figure;
fig.Position(3:4) = [962 673];
k = [0 3 6];
t0 = Time(1);
tf = Time(end);
ax = cell(3,3);
for i=1:3
    k = k + 1;
    ax{i,1} = subplot(3,3,k(1));
    hold on;
    plot(Time, Pref_data(i,:), 'LineWidth',2.0 , 'Color','green', 'LineStyle','-');
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'Color','magenta', 'LineStyle',':');
    plot(Time, P2_data(i,:), 'LineWidth',2.0, 'Color','blue', 'LineStyle','--');
    if (i==1)
        ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
        legend({'$p_{ref}$','$p$','$p_2$'}, 'interpreter','latex', 'fontsize',15, ...
            'Position',[0.3388 0.9427 0.3068 0.0423], 'Orientation','horizontal');
    end

    axis tight;
    hold off;

    ax{i,2} = subplot(3,3,k(2));
    hold on;
    plot(Time, dPref_data(i,:), 'LineWidth',2.0 , 'Color','green', 'LineStyle','-');
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','magenta', 'LineStyle',':');
    plot(Time, dP2_data(i,:), 'LineWidth',2.0, 'Color','blue', 'LineStyle','--');
    for j=1:2, plot([t0 tf], [v_lim(i,j) v_lim(i,j)], 'LineWidth',2, 'Color','red', 'LineStyle','--'); end
    if (i==1), ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15); end

    axis tight;
    hold off;

    ax{i,3} = subplot(3,3,k(3));
    hold on;
    plot(Time, ddPref_data(i,:), 'LineWidth',2.0 , 'Color','green', 'LineStyle','-');
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','magenta', 'LineStyle',':');
    plot(Time, ddP2_data(i,:), 'LineWidth',2.0, 'Color','blue', 'LineStyle','--');
    for j=1:2, plot([t0 tf], [a_lim(i,j) a_lim(i,j)], 'LineWidth',2, 'Color','red', 'LineStyle','--'); end
    if (i==1), ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15); end
    axis tight;
    hold off;
end

%% ===========================================================

function [Timed, yd_data, dyd_data, ddyd_data] = createData(t0, tf, Ts, p0, pf, pm)

    Timed = t0:Ts:tf;
    
    yd_data = get5thOrderPol(p0, pf, Timed);

    tm = [0.6 0.5 0.4]*Timed(end);
    sigma_m = [0.9 1.1 1.2];

    ym = zeros(size(yd_data));
    for i=1:size(ym,1)
       ym(i,:) = pm(i)*exp(-0.5*(Timed-tm(i)).^2/sigma_m(i)^2);
    end
%     ym = pm*exp(-0.5*(Timed-tm).^2/0.9^2);
    
    yd_data = yd_data + ym;
    dyd_data = zeros(size(yd_data));
    ddyd_data = zeros(size(yd_data));
    for i=1:3
        dyd_data(i,:) = [diff(yd_data(i,:)) 0]/Ts; 
        ddyd_data(i,:) = [diff(dyd_data(i,:)) 0]/Ts; 
    end

end

function [P_data, dP_data, ddP_data] = getOutput(gmp, x_data, x_dot, x_ddot)
    
    if (nargin < 4), x_ddot=0; end
    
    n_data = length(x_data);
    n_dof = gmp.numOfDoFs();
    
    P_data = zeros(n_dof, n_data);
    dP_data = zeros(n_dof, n_data);
    ddP_data = zeros(n_dof, n_data);
    
    for j=1:n_data
        x = x_data(j);
        P_data(:,j) = gmp.getYd(x);
        dP_data(:,j) = gmp.getYdDot(x,x_dot);
        ddP_data(:,j) = gmp.getYdDDot(x,x_dot,x_ddot);
    end
    
end
