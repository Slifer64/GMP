clc;
close all;
clear;

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load data
fid = FileIO('data/results.bin', FileIO.in);
data{1} = readData(fid, '1_');
data{2} = readData(fid, '2_');
data{3} = readData(fid, '3_');
data{4} = readData(fid, '4_');
% data{5} = readData(fid, '5_');

pos_lim = fid.read('pos_lim');
vel_lim = fid.read('vel_lim');
accel_lim = fid.read('accel_lim');
y0 = fid.read('y0');
yg = fid.read('yg');
ygd = fid.read('ygd');
tau = fid.read('tau');
view_ = fid.read('view');
  
fid.close();


%% ======== Plot Results ==========

% plot 3D path
fig = figure;
fig.Position(3:4) = [815 716];
ax = axes();
hold on;
plot3(y0(1), y0(2), y0(3), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',12);
plot3(yg(1), yg(2), yg(3), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',12);
plot3(ygd(1), ygd(2), ygd(3), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',12);
legend_ = {};
for k=1:length(data)
    if (~data{k}.plot3D), continue; end
    plot3(data{k}.Pos(1,:), data{k}.Pos(2,:), data{k}.Pos(3,:), 'LineWidth', 3, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    legend_ = [legend_ data{k}.legend];
end
plot3([ygd(1) yg(1)], [ygd(2) yg(2)], [ygd(3) yg(3)], 'LineWidth', 1, 'LineStyle','--', 'Color',[1 0 1 0.5]);
legend(['$p_0$','$g$','$g_d$' legend_], 'interpreter','latex', 'fontsize',17, 'Position',[0.8174 0.6641 0.1531 0.3144]);
axis tight;
x_lim = ax.XLim + 0.05*[-1 1];
y_lim = ax.YLim + 0.05*[-1 1];
z_lim = ax.ZLim + 0.05*[-1 1];
plot3Dbounds(ax, pos_lim)
view(view_);
grid on;
ax.XLim=x_lim; ax.YLim=y_lim; ax.ZLim=z_lim;
ax.FontSize = 14;
hold off;

title_ = {'x coordinate', 'y coordinate', 'z coordinate'};
label_font = 17;
ax_fontsize = 14;

n_dof = length(y0);

for i=1:n_dof
    fig = figure;
    fig.Position(3:4) = [842 1110];

    ax = subplot(3,1,1);
    hold on;
    % plot position trajectory
    legend_ = {};
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
        legend_ = [legend_ data{k}.legend];
    end
    axis tight;
    % plot start and final positions
    plot(0, y0(i), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',10);
    plot(tau, yg(i), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',10);
    plot(tau, ygd(i), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',10); 
    % plot bounds
    plot(ax.XLim, [pos_lim(i,1) pos_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [pos_lim(i,2) pos_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    % labels, title ...
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',label_font);
%     title(title_{i}, 'interpreter','latex', 'fontsize',18);
    legend(legend_, 'interpreter','latex', 'fontsize',17, 'Position',[0.2330 0.9345 0.5520 0.0294], 'Orientation', 'horizontal');
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,2);
    hold on;
    for k=1:length(data) 
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [vel_lim(1) vel_lim(1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [vel_lim(2) vel_lim(2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,3);
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [accel_lim(1) accel_lim(1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [accel_lim(2) accel_lim(2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

end

%% ============================================

function dat = readData(fid, prefix)

    Time = fid.read([prefix 'Time']);
    P_data = fid.read([prefix 'Pos']);
    dP_data = fid.read([prefix 'Vel']);
    ddP_data = fid.read([prefix 'Accel']);
    linestyle = fid.read([prefix 'linestyle']);
    color = fid.read([prefix 'color']);
    legend_ = fid.read([prefix 'legend']);
    plot3D = fid.read([prefix 'plot3D']);
    plot2D = fid.read([prefix 'plot2D']);
    
    dat = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',linestyle, ...
    'color',color, 'legend',legend_, 'plot3D',plot3D, 'plot2D',plot2D);

end

function plot3Dbounds(ax, bounds)
    
    x1 = bounds(1,1);
    x2 = bounds(1,2);
    y1 = bounds(2,1);
    y2 = bounds(2,2);
    z1 = bounds(3,1);
    z2 = bounds(3,2);

%     patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z1 z1 z1 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z2 z2 z2 z2], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

%     patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

end


