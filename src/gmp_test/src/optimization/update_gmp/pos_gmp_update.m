clc;
close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load GMP
gmp = GMP();
gmp_.read(gmp, 'data/gmp_pos.bin','up_');

n_dof = gmp.numOfDoFs();

if (n_dof ~= 3), error('The number of DoFs must be equal to 3!\n'); end

train_err = gmp.autoRetrain(15, 1.5)

%% =============  GMP update  =============
kt = 0.5;
ks = 2;

T = 10/kt;
x_dot = 1/T;
x_ddot = 0;

p0d = gmp.getYd(0); %Pd_data(:,1);
pgd = gmp.getYd(1); %Pd_data(:,end);
p0 = p0d;
pg = ks*(pgd-p0d) + p0;
gmp.setY0(p0);
gmp.setGoal(pg);

t1 = 1.5 / kt;
x1 = t1/T;
p1 = gmp.getYd(x1) + [-0.1; 0.08; -0.12];

t2 = 3 / kt;
x2 = t2/T;
p2_dot = gmp.getYdDot(x2,x_dot) + [0.05; -0.08; 0.1];

t3 = 4.25 / kt;
x3 = t3/T;
p3_ddot = gmp.getYdDDot(x3,x_dot,x_ddot) + [-0.02; 0.03; 0.05];

t4 = 5.5 / kt;
x4 = t4/T;
p4 = gmp.getYd(x4) + [0.2; -0.15; -0.13].*[1; 1; 1];
p4_ddot = gmp.getYdDDot(x4,x_dot,x_ddot) + [0.02; -0.04; -0.025];
s4 = [GMP_phase(x4,x_dot,x_ddot), GMP_phase(x4,x_dot,x_ddot)];
type4 = [GMP_UpdateType.POS, GMP_UpdateType.ACCEL];
Z4 = [p4, p4_ddot];

t5 = 7 / kt;
x5 = t5/T;
p5_dot = gmp.getYdDot(x5,x_dot) + [0.1; 0.05; -0.1].*[1; 1; 1];
p5_ddot = gmp.getYdDDot(x5,x_dot,x_ddot) + [-0.03; 0.02; 0.025];
s5 = [GMP_phase(x5,x_dot,x_ddot), GMP_phase(x5,x_dot,x_ddot)];
type5 = [GMP_UpdateType.VEL, GMP_UpdateType.ACCEL];
Z5 = [p5_dot, p5_ddot];

% make a deep copy of the original GMP
gmp0 = gmp.deepCopy();

gmp_up = GMP_Update(gmp);
gmp_up.initSigmaWfromMsr(0:0.01:1);
gmp_up.enableSigmawUpdate(true);
gmp_up.setMsrNoiseVar(1e-4);

O_vec = zeros(3,1);

% make sure boundary conditions are satisfied
gmp_up.updatePos(0, p0);
gmp_up.updatePos(1, pg);
gmp_up.updateVel(0, x_dot, O_vec);
gmp_up.updateVel(1, x_dot, O_vec);
gmp_up.updateAccel(0, x_dot, x_ddot, O_vec);
gmp_up.updateAccel(1, x_dot, x_ddot, O_vec);

gmp_up.updatePos(x1, p1);
gmp_up.updateVel(x2, x_dot, p2_dot);
gmp_up.updateAccel(x3, x_dot, x_ddot, p3_ddot);
gmp_up.updateWeights(s4, Z4, type4);
gmp_up.updateWeights(s5, Z5, type5);


%% Generate original and new trajectories

Ts = 0.01;
Time = (0:Ts:T);
x = Time/Time(end);
 
N = length(x);
P_new_data = zeros(3,N);
dP_new_data = zeros(3,N);
ddP_new_data = zeros(3,N);

P_data = zeros(3,N);
dP_data = zeros(3,N);
ddP_data = zeros(3,N);

for j=1:N
    P_data(:,j) = gmp0.getYd(x(j));
    dP_data(:,j) = gmp0.getYdDot(x(j), x_dot);
    ddP_data(:,j) = gmp0.getYdDDot(x(j), x_dot, x_ddot);
    
    P_new_data(:,j) = gmp.getYd(x(j));
    dP_new_data(:,j) = gmp.getYdDot(x(j), x_dot);
    ddP_new_data(:,j) = gmp.getYdDDot(x(j), x_dot, x_ddot);
end


%% Plot results

points(1) = getPoint(t1, x1, x_dot, x_ddot, p1, [], []);
points(2) = getPoint(t2, x2, x_dot, x_ddot, [], p2_dot, []);
points(3) = getPoint(t3, x3, x_dot, x_ddot, [], [], p3_ddot);
points(4) = getPoint(t4, x4, x_dot, x_ddot, p4, [], p4_ddot);
points(5) = getPoint(t5, x5, x_dot, x_ddot, [], p5_dot, p5_ddot);

fig = figure;
fig.Position(3:4) = [980 756];
k = [0 3 6];
ax = cell(3,3);
title_ = {'$x$', '$y$', '$z$'};
for i=1:3
    k = k + 1;
    ax{i,1} = subplot(3,3,k(1));
    hold on;
    plot(Time, P_new_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ax{i,1}.FontSize = 14;
    title(title_{i}, 'interpreter','latex', 'fontsize',19);
    if (i==1), ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',19); end
    if (i==1)
        legend({'modified','original'}, 'interpreter','latex', 'fontsize',18, ...
        'Position',[0.4008 0.9623 0.2528 0.0440], 'Orientation','horizontal');
    end
    axis tight;
    hold off;

    ax{i,2} = subplot(3,3,k(2));
    hold on;
    plot(Time, dP_new_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ax{i,2}.FontSize = 14;
    if (i==1), ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',19); end
    axis tight;
    hold off;

    ax{i,3} = subplot(3,3,k(3));
    hold on;
    plot(Time, ddP_new_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    if (i==1), ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',19); end
    ax{i,3}.FontSize = 14;
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',18);
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
