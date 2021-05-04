clc;
% close all;
clear;

addpath('../../matlab/lib/io_lib/');
import_io_lib();

addpath('../../matlab/lib/gmp_lib/');
import_gmp_lib();

%% Load training data
gmp = GMP_nDoF(1, 2);
gmp_new = GMP_nDoF(1, 2);

GMP_nDoF_IO.read(gmp, 'data/gmp_pos.bin','up_');
GMP_nDoF_IO.read(gmp_new, 'data/gmp_pos_updated.bin','up_');

kt = 0.5;
ks = 2;
p0d = gmp.getYd(0); %Pd_data(:,1);
pgd = gmp.getYd(1); %Pd_data(:,end);
p0 = p0d;
pg = ks*(pgd-p0d) + p0;
T = 10/kt;
x_dot = 1/T;
x_ddot = 0;

gmp.setY0(p0);
gmp.setGoal(pg);

gmp_new.setY0(p0);
gmp_new.setGoal(pg);


points = getPoint([], [], [], [], [], [], []);

t1 = 1.5 / kt;
x1 = t1/T;
p1 = gmp.getYd(x1) + [-0.1; 0.08; -0.12].*[1; 1; 1];
points(1) = getPoint(t1, x1, x_dot, x_ddot, p1, [], []);

t2 = 3 / kt;
x2 = t2/T;
p2_dot = gmp.getYdDot(x2,x_dot) + [0.05; -0.08; 0.1].*[1; 1; 1];
points(2) = getPoint(t2, x2, x_dot, x_ddot, [], p2_dot, []);

t3 = 4.25 / kt;
x3 = t3/T;
p3_ddot = gmp.getYdDDot(x3,x_dot,x_ddot) + [-0.02; 0.03; 0.05].*[1; 1; 1];
points(3) = getPoint(t3, x3, x_dot, x_ddot, [], [], p3_ddot);

t4 = 5.5 / kt;
x4 = t4/T;
p4 = gmp.getYd(x4) + [0.2; -0.15; -0.13].*[1; 1; 1];
p4_ddot = gmp.getYdDDot(x4,x_dot,x_ddot) + [0.02; -0.04; -0.025].*[1; 1; 1];
points(4) = getPoint(t4, x4, x_dot, x_ddot, p4, [], p4_ddot);

t5 = 7 / kt;
x5 = t5/T;
p5_dot = gmp.getYdDot(x5,x_dot) + [0.1; 0.05; -0.1].*[1; 1; 1];
p5_ddot = gmp.getYdDDot(x5,x_dot,x_ddot) + [-0.03; 0.02; 0.025].*[1; 1; 1];
points(5) = getPoint(t5, x5, x_dot, x_ddot, [], p5_dot, p5_ddot);


%% simulate
Ts = 0.01;
Time = 0:Ts:T;
x = Time/Time(end);
n_data = length(x);

P_data = zeros(3,n_data);
dP_data = zeros(3,n_data);
ddP_data = zeros(3,n_data);

P_new_data = zeros(3,n_data);
dP_new_data = zeros(3,n_data);
ddP_new_data = zeros(3,n_data);
for j=1:n_data
    P_data(:,j) = gmp.getYd(x(j));
    dP_data(:,j) = gmp.getYdDot(x(j), x_dot);
    ddP_data(:,j) = gmp.getYdDDot(x(j), x_dot, x_ddot);
    
    P_new_data(:,j) = gmp_new.getYd(x(j));
    dP_new_data(:,j) = gmp_new.getYdDot(x(j), x_dot);
    ddP_new_data(:,j) = gmp_new.getYdDDot(x(j), x_dot, x_ddot);
end


%% Plot results
fig = figure;
k = [0 3 6];
ax = cell(3,3);
for i=1:3
    k = k + 1;
    ax{i,1} = subplot(3,3,k(1));
    hold on;
    plot(Time, P_new_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    legend({'modified','original'}, 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    ax{i,2} = subplot(3,3,k(2));
    hold on;
    plot(Time, dP_new_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    ax{i,3} = subplot(3,3,k(3));
    hold on;
    plot(Time, ddP_new_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end

for i=1:length(points), plotUpdatePoint(points(i), ax); end

%% ==================================================================
%% ==================================================================

function point = getPoint(t, x, x_dot, x_ddot, p, p_dot, p_ddot)

    point = struct('t',t, 'x',x, 'x_dot',x_dot, 'x_ddot',x_ddot, 'p',p, 'p_dot',p_dot, 'p_ddot',p_ddot);

end

function plotUpdatePoint(s, ax)

    if (isempty(s)), return; end

    for i=1:3

        if (~isempty(s.p))
            ax_ = ax{i,1};
            hold(ax_, 'on');
            scatter([s.t], [s.p(i)], 'Marker','*', 'MarkerEdgeColor',[1 0 0], 'MarkerEdgeAlpha',0.7, 'LineWidth',2, 'SizeData', 100, 'Parent',ax_, 'HandleVisibility','off');
            plot([s.t s.t], ax_.YLim, 'LineWidth',1, 'Color','green', 'LineStyle','--', 'Parent',ax_, 'HandleVisibility','off');
            hold(ax_, 'off');
        end

        if (~isempty(s.p_dot))
            ax_ = ax{i,2};
            hold(ax_, 'on');
            scatter([s.t], [s.p_dot(i)], 'Marker','*', 'MarkerEdgeColor',[1 0 0], 'MarkerEdgeAlpha',0.7, 'LineWidth',2, 'SizeData', 100, 'Parent',ax_, 'HandleVisibility','off');
            plot([s.t s.t], ax_.YLim, 'LineWidth',1, 'Color','green', 'LineStyle','--', 'Parent',ax_, 'HandleVisibility','off');
            hold(ax_, 'off');
        end

        if (~isempty(s.p_ddot))
            ax_ = ax{i,3};
            hold(ax_, 'on');
            scatter([s.t], [s.p_ddot(i)], 'Marker','*', 'MarkerEdgeColor',[1 0 0], 'MarkerEdgeAlpha',0.7, 'LineWidth',2, 'SizeData', 100, 'Parent',ax_, 'HandleVisibility','off');
            plot([s.t s.t], ax_.YLim, 'LineWidth',1, 'Color','green', 'LineStyle','--', 'Parent',ax_, 'HandleVisibility','off');
            hold(ax_, 'off');
        end

    end
end
