clc;
% close all;
clear;

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

%% Load data
fid = FileIO('pos_results.bin', FileIO.in);
Time         = fid.read('Time');
P_data       = fid.read('P_data');
dP_data      = fid.read('dP_data');
ddP_data     = fid.read('ddP_data');
P_new_data   = fid.read('P_new_data');
dP_new_data  = fid.read('dP_new_data');
ddP_new_data = fid.read('ddP_new_data');

points = getPoint([], [], [], [], [], [], []);

gmp = GMP();
gmp_.read(gmp, 'gmp_pos.bin','up_');

kt = 0.5;
T = 10/kt;
x_dot = 1/T;
x_ddot = 0;
  
t1 = 1.5 / kt;
x1 = t1/T;
p1 = fid.read('p1');

t2 = 3 / kt;
x2 = t2/T;
p2_dot = fid.read('p2_dot');

t3 = 4.25 / kt;
x3 = t3/T;
p3_ddot = fid.read('p3_ddot');

t4 = 5.5 / kt;
x4 = t4/T;
p4 = fid.read('p4');
p4_ddot = fid.read('p4_ddot');

t5 = 7 / kt;
x5 = t5/T;
p5_dot = fid.read('p5_dot');
p5_ddot = fid.read('p5_ddot');


%% Plot results

points(1) = getPoint(t1, x1, x_dot, x_ddot, p1, [], []);
points(2) = getPoint(t2, x2, x_dot, x_ddot, [], p2_dot, []);
points(3) = getPoint(t3, x3, x_dot, x_ddot, [], [], p3_ddot);
points(4) = getPoint(t4, x4, x_dot, x_ddot, p4, [], p4_ddot);
points(5) = getPoint(t5, x5, x_dot, x_ddot, [], p5_dot, p5_ddot);

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
