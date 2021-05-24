clc;
% close all;
clear;

set_matlab_utils_path('../');

%% Load training data
% load('data/train_data.mat', 'Data');
%
% Timed = Data.Time;
% Pd_data = Data.Pos;
% dPd_data = Data.Vel;
% ddPd_data = Data.Accel;

Timed = 0:0.002:10;
[Pd_data, dPd_data, ddPd_data] = get5thOrderPol([0 0 0]', [0.8 0.7 0.9]', Timed);

Ts = Timed(2)-Timed(1);

%% initialize and train GMP
train_method = 'LS';
N_kernels = 20;
kernels_std_scaling = 1;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

% x = Timed/Timed(end);
% gmp.plotPsi(x);

% gmp.plotWeightsCovariance();

% GMP_IO.write(gmp, 'gmp_pos.bin','up_');
% return

%% scale spatio-temporally
kt = 0.5;
ks = 2;
p0d = gmp.getYd(0); %Pd_data(:,1);
pgd = gmp.getYd(1); %Pd_data(:,end);
p0 = p0d;
pg = ks*(pgd-p0d) + p0;
T = Timed(end)/kt;
Time = Timed/kt;
x = Time / T;
x_dot = 1/T;
x_ddot = 0;

gmp.setY0(p0);
gmp.setGoal(pg);

N = length(x);
P_data = zeros(n_dof,N);
dP_data = zeros(n_dof,N);
ddP_data = zeros(n_dof,N);

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

gmp_up = GMP_Update(gmp);
% gmp_up.initSigmaWfromMsr(0:0.01:1);
gmp_up.enableSigmawUpdate(true);
gmp_up.setMsrNoiseVar(1e-4);
for i=1:length(points), updateGMP(gmp_up, points(i)); end
% updateGMP_all(gmp_up, points);


%% simulate
for j=1:N
    P_data(:,j) = gmp.getYd(x(j));
    dP_data(:,j) = gmp.getYdDot(x(j), x_dot);
    ddP_data(:,j) = gmp.getYdDDot(x(j), x_dot, x_ddot);
end

% for i=1:3
%     dP_data(i,:) = [diff(P_data(i,:)) 0] / Ts * kt;
%     ddP_data(i,:) = [diff(dP_data(i,:)) 0] / Ts * kt;
% end


%% obtain scaled demo data
Timed = Timed / kt;
Pd_data = ks*(Pd_data-p0d) + p0;
dPd_data = ks*kt*dPd_data;
ddPd_data = ks*kt^2*ddPd_data;


%% Plot results
fig = figure;
k = [0 3 6];
ax = cell(3,3);
for i=1:3
    k = k + 1;
    ax{i,1} = subplot(3,3,k(1));
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Timed, Pd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    legend({'modified','original'}, 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    ax{i,2} = subplot(3,3,k(2));
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, dPd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    ax{i,3} = subplot(3,3,k(3));
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, ddPd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
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

function updateGMP(gmp_up, point)

    updateGMP_all(gmp_up, point);

end

function updateGMP_all(gmp_up, points)

    if (isempty(points)), return; end

    s = [];
    type = [];
    z = [];

    for k=1:length(points)
        point = points(k);
        si = GMP_phase(point.x, point.x_dot, point.x_ddot);

        if (~isempty(point.p))
            s = [s si];
            z = [z point.p];
            type = [type 0];
        end

        if (~isempty(point.p_dot))
            s = [s si];
            z = [z point.p_dot];
            type = [type 1];
        end

        if (~isempty(point.p_ddot))
            s = [s si];
            z = [z point.p_ddot];
            type = [type 2];
        end
    end

    gmp_up.updateWeights(s, z, type);

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
