clc;
close all;
clear;

%% =============  includes...  =============
import_gmp_lib();
import_io_lib();

%% =============  Create training data  =============

% Timed = 0:0.005:6;
% [Yd_data, Yd_dot_data, Yd_ddot_data] = get5thOrderPol([0 0 0]', [0.3 0.5 0.8]', Timed);

load('train_data.mat','data');
Timed = data.Time;
Yd_data = data.Pos;


y0 = Yd_data(:,1);
yg0 = Yd_data(:,end);

n_dof = size(Yd_data, 1);


%% =============  Create/Train GMP  =============

train_method = 'LS';
N_kernels = 25;
kernels_std_scaling = 1.5;

gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Yd_data);
offline_train_mse
toc


%% =============  Simulate GMP  =============
y = y0;
y_dot = zeros(n_dof,1);
y_ddot = zeros(n_dof,1);

Tf = Timed(end);

dt = 0.002;

t = 0;
s = 0;
s_dot = 1/Tf;
s_ddot = 0;

Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];

ks = diag([1.4 0.8 1]);
yg = y0 + ks*(yg0-y0);

gmp.setY0(y0);
gmp.setGoal(yg);

gmp_up = GMP_Update(gmp);
% i1 = 3;
% i2 = 8;
% i3 = i1 + (i2-i1)/2;
% for i=i1:0.5:i2
%     s1 = i/10;
%     y1 = gmp.getYd(s1) + [0.03*(i2 - abs(i-i3)); 0; 0];
%     gmp_up.updatePos(s1, y1);
% end

obst = struct('c',[0.49; 0.26; 0.51], 'r',0.05);

while (true)

    if (s >= 1), break; end

    yd = gmp.getYd(s);
    yd_dot = gmp.getYdDot(s,s_dot);
    yd_ddot = gmp.getYdDDot(s,s_dot,s_ddot);

    y_ddot = -50*(y_dot - yd_dot) - 300*(y - yd) + yd_ddot;

    if ( norm(y-obst.c) <= 1*obst.r )
        % find direction of motion
        v1 = yd_dot/norm(yd_dot);
        v2 = cross([0; 0; 1], v1);

        s1 = s;
        y1 = y + 0.2*v2 + obst.r*v1;
        s2 = s1 + 0.1;
        y2 = y - 0.2*v2 + obst.r*v1;

        gmp_up.updatePos(s1, y1);
        gmp_up.updatePos(s2, y2);

    end

    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data y_dot];
    ddY_data = [ddY_data y_ddot];

    t = t + dt;
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;
    y = y + y_dot*dt;
    y_dot = y_dot + y_ddot*dt;

end

% data = struct('Time',Time, 'Pos',Y_data, 'Vel',dY_data, 'Accel',ddY_data);
% save('train_data.mat','data');

% vel_norm = arrayfun( @(j) norm(dY_data(:,j)) , 1:size(dY_data,2), 'uni','false');
vel_norm = zeros(1, length(Time));
for j=1:size(dY_data,2), vel_norm(j) = norm(dY_data(:,j)); end
% figure;
% plot(Time, vel_norm, 'LineWidth',2, 'Color','magenta');

fig = figure;
fig.Position(3:4) = [944 666];
ax = axes('NextPlot','add');
plot3(Y_data(1,:), Y_data(2,:), Y_data(3,:), 'LineWidth',2);
plot3(Yd_data(1,:), Yd_data(2,:), Yd_data(3,:), 'LineWidth',2, 'LineStyle','--');
plot3(y0(1), y0(2), y0(3), 'LineStyle','none', 'Marker','o', 'MarkerSize',10, 'LineWidth',4, 'Color',[0 0.8 0]);
plot3(yg0(1), yg0(2), yg0(3), 'LineStyle','none', 'Marker','x', 'MarkerSize',10, 'LineWidth',4, 'Color','magenta');
plot3(yg(1), yg(2), yg(3), 'LineStyle','none', 'Marker','x', 'MarkerSize',10, 'LineWidth',4, 'Color',[0.8 0 0]);

plotSphere(ax, obst.r, obst.c);

xlabel('$X$', 'interpreter','latex', 'fontsize',15);
ylabel('$Y$', 'interpreter','latex', 'fontsize',15);
zlabel('$Z$', 'interpreter','latex', 'fontsize',15);
grid on;
view(-19.8, 16.3);


%% ===================================================================

function s = plotSphere(ax, r, c)

    [X, Y, Z] = sphere(10);

    X = r*X + c(1);
    Y = r*Y + c(2);
    Z = r*Z + c(3);

    s = surf(X,Y,Z, 'Parent',ax, 'EdgeColor','none', 'FaceColor',[0.85 0.33 0.1], 'FaceAlpha',0.5);

end

function [y, y_dot, y_ddot] = get5thOrderPol(y0, yf, Time)

    y0 = y0(:);
    yf = yf(:);
    Time = Time(:)';

    T = Time(end);
    t = Time/T;

    y = y0 + (yf - y0) * (10*t.^3 - 15*t.^4 + 6*t.^5 );

    if (nargout > 1)
        y_dot = (yf - y0) * (30*t.^2 - 60*t.^3 + 30*t.^4 ) / T;
    end

    if (nargout > 2)
        y_ddot = (yf - y0) * (60*t - 180*t.^2 + 120*t.^3 ) / T^2;
    end

end
