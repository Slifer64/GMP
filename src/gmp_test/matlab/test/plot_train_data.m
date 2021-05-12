clc;
close all;
clear;

set_matlab_utils_path();

data = loadData('data/train_data.bin');
data2 = loadData('data/train_data2.bin');

n_data = length(data2.Time);

[x_data, fv_data] = loadPhaseFvData('data/x_fv_data.bin', n_data);

ticks_fs = 13;
axis_color = {[0 0.45 0.74]; [0.85 0.33 0.1]; [0.93 0.69 0.13]};
legend_fs = 17;
axis_lb_fs = 17; 
title_fs = 18;

%% ====================================================
fig = figure;
fig.Position(3:4) = [524 723];
% --------------------------------------------
ax = subplot(4,1,1); hold on; ax.FontSize = ticks_fs;
for i=1:3
%     plot(data2.Time, data2.Pos(i,:), 'LineWidth',2.5, 'Color',axis_color{i});
%     plot(data.Time, data.Pos(i,:), 'LineWidth',2.5, 'Color',axis_color{i}, 'LineStyle','--');
    plot(data.Time, data.Vel(i,:), 'LineWidth',2.5, 'Color',axis_color{i}, 'LineStyle','-');
    ylabel('$m/s$', 'interpreter','latex', 'fontsize',axis_lb_fs);
end
axis tight;
legend({'$\dot{x}$','$\dot{y}$','$\dot{z}$'}, 'interpreter','latex', 'fontsize',legend_fs, 'Position',[0.8487 0.6847 0.1030 0.1198]);
title('Phase 1: Velocity profile', 'interpreter','latex', 'fontsize',title_fs);
% --------------------------------------------
ax = subplot(4,1,2); hold on; ax.FontSize = ticks_fs;
for i=1:3
    plot(data2.Time, data2.Vel(i,:), 'LineWidth',2.5, 'Color',axis_color{i});
end
ylabel('$m/s$', 'interpreter','latex', 'fontsize',axis_lb_fs);
axis tight;
title('Phase 2: Velocity profile', 'interpreter','latex', 'fontsize',title_fs);
% --------------------------------------------
ax = subplot(4,1,3); hold on; ax.FontSize = ticks_fs;
plot(data2.Time, x_data, 'LineWidth',2.5, 'Color','blue', 'DisplayName','x: phase 2');
plot([data2.Time(1) data2.Time(end)], [0 data2.Time(end)/data.Time(end)], 'LineWidth',2, 'Color','magenta', 'LineStyle',':', 'DisplayName','x: phase 1');
legend({}, 'interpreter','latex', 'fontsize',legend_fs, 'Position',[0.1521 0.3992 0.2743 0.0815]);
ylabel('phase var.', 'interpreter','latex', 'fontsize',axis_lb_fs);
axis tight;
% --------------------------------------------
ax = subplot(4,1,4); ax.FontSize = ticks_fs;
plot(data2.Time, fv_data, 'LineWidth',2.5, 'Color','red');
legend({'$f_v$'}, 'interpreter','latex', 'fontsize',legend_fs, 'Position',[0.7675 0.2081 0.1173 0.0432]);
xlabel('time $s$', 'interpreter','latex', 'fontsize',axis_lb_fs);
axis tight;


%% ====================================================
figure;
ax = plot_3Dpath_with_orientFrames(data2.Pos, data2.Quat, 'numberOfFrames',6, 'LineWidth',3, 'frameLineWidth',2, 'frameScale',0.3);
ax.FontSize = 12;
grid on;

%% =============================================================

           
function data = loadData(path)

    fid = FileIO(path, FileIO.in);
    
    Time = fid.read('Timed');
    Pos = fid.read('Pd_data');
    Vel = fid.read('dPd_data');
    Quat = fid.read('Qd_data');
    RotVel = fid.read('vRotd_data');
    JointPos = fid.read('joint_pos_data');
    
    data = struct('Time',Time, 'Pos',Pos, 'Vel',Vel, 'Quat',Quat, 'RotVel',RotVel, 'JointPos',JointPos);

end


function [x_data, fv_data] = loadPhaseFvData(path, len)

    fid = FileIO(path, FileIO.in);
    x_data = fid.read('x_data');
    fv_data = fid.read('fv_data');
    
    n_data = length(x_data);
    
    j1 = 0;
    for j=1:n_data
        if (x_data(j) > 0)
            j1 = j;
            break
        end
    end
    
    j2 = n_data;
    for j=n_data:-1:j1
        if (x_data(j) < 1)
            j2 = j;
            break
        end
    end
    
    x_data = x_data(j1:j2);
    fv_data = fv_data(j1:j2);
    
    n_pad = len - length(x_data);
    if (n_pad > 0)
       x_data = [x_data ones(1,n_pad)*x_data(end)];
       fv_data = [fv_data ones(1,n_pad)*fv_data(end)];
    else
        x_data = x_data(1:len);
        fv_data = fv_data(1:len);
    end

end

    
    