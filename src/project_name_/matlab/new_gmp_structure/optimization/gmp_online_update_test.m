clc;
close all;
clear;

rng(0);

set_matlab_utils_path('../');

%% Generate data
Ts = 0.002;
Timed = 0:Ts:10;
% nominal data
[yd_data, dyd_data, ddyd_data] = get5thOrderPol(0, 0.8, Timed);
y_demo_data = yd_data;
x_demo_data = Timed / Timed(end);

% % pad extra data
% t_end = 16;
% T_extra = Time(end)+Ts:Ts:t_end;
% n_pad = length(T_extra);
% Time = [Time T_extra];
% yd_data = [yd_data repmat(yd_data(end), 1, n_pad)];
% dyd_data = [dyd_data repmat(dyd_data(end), 1, n_pad)];
% ddyd_data = [ddyd_data repmat(ddyd_data(end), 1, n_pad)];


n_data = length(Timed);
tm = Timed(end)/2;
% modified data
ym = 0.25*exp(-0.5*(Timed-tm).^2/0.9^2);
% for j=1:n_data
%    if (norm(ym(:,j)) < 1e-10), ym(:,j) = 0*ym(:,j); end
% end
yh_data = yd_data + ym;
dyh_data = [diff(yh_data) 0]/Ts; 
ddyh_data = [diff(dyh_data) 0]/Ts; 
    
tau = Timed(end);
x_data = Timed / tau;
x_dot = 1/tau;
x_ddot = 0;


%% Initialize GMP
N_kernels = 25;
n_dofs = size(yd_data,1);
gmp = GMP_nDoF(n_dofs, N_kernels, 2);

%% find the points where Pd and Ph differ
ind2 = findUpdatePoints(yd_data, yh_data, 1e-5);
y2 = yh_data(:,ind2);
Time2 = Timed(ind2);
x2 = x_data(ind2);


%% update Weights with LS
[Wo, Swo] = batch_LS(gmp, x_data, yd_data);
% [y0_data, dy0_data, ddy0_data] = getOutput(gmp, x_data, x_dot);


%% update Weights with RLS
% W0 = zeros(n_dofs, N_kernels);
% Sw0 = 1e6*eye(N_kernels,N_kernels);
% R = 1;
% [W, Sw, eig_Sw] = recursive_LS(gmp, x_data, yd_data, W0, Sw0, R);
% [y_data, dy_data, ddy_data] = getOutput(gmp, x_data, x_dot);
% plotEigSw(eig_Sw, [0.95 0.05]);

% w_err = norm(W-Wo)
% err_Sw = norm(Sw(:) - Swo(:))
% plotSw(Sw);
% plotSw(Swo);

%% update Weights with recursive U/D LS
% gmp_temp = gmp.deepCopy();
% y_old = yd_data(ind2);
% % [W2, Sw2, eig_Sw2] = recursive_UD_LS(y2, y_old, H2, W, Sw, 1e-3, 1);
% % [W2, Sw2, eig_Sw2] = recursive_LS(y2, H2, W, Sw, 1e-3);
% [W2, Sw2, eig_Sw2] = recursive_LS(gmp_temp, x2, y2, W, Sw, 1e-3);
% [y2_data, dy2_data, ddy2_data] = getOutput(gmp_temp, x_data, x_dot);
% % plotEigSw(eig_Sw2, [0.95 0.95]);


%% on-line update with external force

Swo = 1e-1*eye(N_kernels, N_kernels);
% plotSw(Swo)
% return 

gmp.W = Wo;
gmp_up = GMP_nDoF_Update(gmp);
gmp_up.setSigmaW(Swo);
gmp_up.setMsrNoiseVar(1);
gmp_up.enableSigmawUpdate(true);

t_end = 15;
dt = Ts;
Time = 0:dt:t_end;


n_data = length(Time);
n0 = length(Timed);

yh_data = [yh_data repmat(yh_data(:,end), 1,n_data-n0)];
dyh_data = [dyh_data repmat(dyh_data(:,end), 1,n_data-n0)];
ddyh_data = [ddyh_data repmat(ddyh_data(:,end), 1,n_data-n0)];
x_data = [x_data repmat(x_data(end), 1,n_data-n0)];


H = zeros(N_kernels, n_data);
H_dot = zeros(N_kernels, n_data);
H_ddot = zeros(N_kernels, n_data);
for j=1:length(x_data)
    x = x_data(j);
    H(:,j) = gmp.regressVec(x);
    H_dot(:,j) = gmp.regressVecDot(x, x_dot);
    H_ddot(:,j) = gmp.regressVecDDot(x, x_dot, x_ddot); 
end


yref_data = []; %zeros(1,n_data);
dyref_data = []; %zeros(1,n_data);
ddyref_data = []; %zeros(1,n_data);
y_data = []; %zeros(1,n_data);
dy_data = []; %zeros(1,n_data);
ddy_data = []; %zeros(1,n_data);
Fext_data = []; %zeros(1,n_data);

y = yd_data(1);
dy = 0;
ddy = 0;

v = zeros(n_dofs,1);
% v_data = zeros(1,n_data);
% for j=1:n_data
%    v_data(j) = v;
%    v = 0.99*v + 0.1*randn(n_dofs,1);
% end
% figure;
% plot(Time, v_data);
% return

mult = [100 1 1];
fig = figure;
fig.Position(3:4) = [697 751];
sh = cell(4,1);
for i=1:3
    ax = subplot(4,1,i); hold(ax,'on');
    pl_y = plot(Time, nan*ones(1,n_data), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
    pl_y2 = plot(Time, nan*ones(1,n_data), 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
    pl_yref = plot(Time, nan*ones(1,n_data), 'LineWidth',2, 'Color','red', 'LineStyle','-.');
    sh{i} = struct('ax',ax, 'pl_y',pl_y, 'pl_y2',pl_y2, 'pl_yref',pl_yref);
    if (i==1)
        plot(Time, mult(i)*yh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
        lg = legend({'$y$','$y_2$','$y_{ref}$','$y_h$'}, 'interpreter','latex', 'fontsize',15, ...
            'Position',[0.2563 0.9385 0.5357 0.0354], 'Orientation','horizontal');
    elseif (i==2)
        plot(Time, mult(i)*dyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
    else
        plot(Time, mult(i)*ddyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
    end
end
subplot(4,1,4);
pl_Fext = plot(Time, zeros(1,n_data), 'LineWidth',2, 'Color','red');
plot_count = 0;

for j=1:n_data
    
    t = Time(j);
    x = x_data(j);
    
    if (x >= 1)
       x = 1;
       x_dot = 0;
       x_ddot = 0;
    end

    % get reference
    y_ref = gmp.getYd(x);
    
    % update weights
    z = y;
    %z = yh_data(:,j);
    if (norm(z - y_ref) > 1e-5)
        
        s = repmat(GMP_phase(x,x_dot,x_ddot), 1,1);
        Z = [y];
        type = [GMP_UpdateType.POS];
        r_n = [1e-3];
        
%         s = repmat(GMP_phase(x,x_dot,x_ddot), 1,2);
%         Z = [y dy];
%         type = [GMP_UpdateType.POS, GMP_UpdateType.VEL];
%         r_n = [1e-3, 1e-2];
        
%         s = repmat(GMP_phase(x,x_dot,x_ddot), 1,3);
%         Z = [y dy ddy];
%         type = [GMP_UpdateType.POS, GMP_UpdateType.VEL, GMP_UpdateType.ACCEL];
%         r_n = [1e-3, 1e-2, 5e-2];
        
        gmp_up.updateWeights(s, Z, type, r_n);

%         gmp_up.updatePos(x, z);
    end
    
    y_ref = gmp.getYd(x);
    dy_ref = gmp.getYdDot(x, x_dot);
    ddy_ref = gmp.getYdDDot(x, x_dot, x_ddot);

    v = 0.99*v + 0.05*randn(n_dofs,1);
    if ((t-Timed(end))<1), v = 0*v; end
 
    Fext = 350*(yh_data(:,j) - y) + 0*80*(dyh_data(:,j) - dy) + 0*ddyh_data(:,j); % + v;
    %Fc = 
    
    ddy = ddy_ref + ( 30*(dy_ref - dy) + 200*(y_ref - y) + Fext);
    
    % logging
    yref_data(:,j) = y_ref;
    dyref_data(:,j) = dy_ref;
    ddyref_data(:,j) = ddy_ref;
    y_data(:,j) = y;
    dy_data(:,j) = dy;
    ddy_data(:,j) = ddy;
    Fext_data(:,j) = Fext;
    
    % Plot on-line
    s_y = [y dy ddy];
    s_y2 = [gmp.W*H; gmp.W*H_dot; gmp.W*H_ddot];
    s_yref = [y_ref dy_ref ddy_ref];
    for i=1:3
        sh{i}.pl_y.YData(j) = mult(i)*s_y(:,i);
        sh{i}.pl_y2.YData = mult(i)*s_y2(i,:);
        sh{i}.pl_yref.YData(j) = mult(i)*s_yref(:,i);
    end
    pl_Fext.YData(j) = abs(Fext);
    plot_count = plot_count + 1;

    if (plot_count == 100)
        plot_count = 0;
        drawnow
        pause
    end
    
    % integration
    y = y + dy*dt;
    dy = dy + ddy*dt;
    
end

[y2_data, dy2_data, ddy2_data] = getOutput(gmp, x_data, x_dot);

figure;
for i=1:n_dofs
   subplot(n_dofs,1,i);
   plot(Time, abs(Fext_data(i,:)), 'LineWidth',2, 'Color','red');
   axis tight;
end


%% plot trajectories
figure;
subplot(3,1,1); hold on;
plot(Time, y_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, y2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
plot(Time, yref_data, 'LineWidth',2, 'Color','red', 'LineStyle','-.');
% plot(Time, yd_data, 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
plot(Time, yh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
% plot(Time2, y2, 'LineWidth',2.0, 'LineStyle',':', 'Color',[0 1 0 0.6]);
% ylim([min(yh_data) max(yh_data)]);
subplot(3,1,2); hold on;
plot(Time, dy_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, dy2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
plot(Time, dyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
plot(Time, dyref_data, 'LineWidth',2, 'Color','red', 'LineStyle','-.');
ylim([min(dyh_data) max(dyh_data)]);
subplot(3,1,3); hold on;
plot(Time, ddy_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, ddy2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
plot(Time, ddyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
plot(Time, ddyref_data, 'LineWidth',2, 'Color','red', 'LineStyle','-.');
ylim([min(ddyh_data) max(ddyh_data)]);

%% ===========================================================

%% Trajectory output functions

function ind = findUpdatePoints(yd_data, yh_data, tol)

    ind = [];
    n_data = length(yd_data);
    for i=1:n_data
       if (norm(yd_data(i)-yh_data(i)) > tol), ind = [ind i]; end
    end
    
end

function [y_data, dy_data, ddy_data] = getOutput(gmp, x_data, x_dot, x_ddot)
    
    if (nargin < 4), x_ddot=0; end
    
    n_data = length(x_data);
    n_dof = gmp.numOfDoFs();
    
    y_data = zeros(n_dof, n_data);
    dy_data = zeros(n_dof, n_data);
    ddy_data = zeros(n_dof, n_data);
    
    for j=1:n_data
        x = x_data(j);
        y_data(:,j) = gmp.getYd(x);
        dy_data(:,j) = gmp.getYdDot(x,x_dot);
        ddy_data(:,j) = gmp.getYdDDot(x,x_dot,x_ddot);
    end
    
end


%% Update weights functions

function [W, Sw] = batch_LS(gmp, x_data, yd_data)

    [~, Sw] = gmp.train('LS', x_data, yd_data);
    W = gmp.W;

end

function [W, Sw, eig_Sw] = recursive_LS(gmp, x, y, W0, Sw0, R)

    if (nargin < 6), R=1; end
    gmp.W = W0;
    
    gmp_up = GMP_nDoF_Update(gmp);
    gmp_up.setSigmaW(Sw0);
    gmp_up.setMsrNoiseVar(R);
    gmp_up.enableSigmawUpdate(true);
    
    n = size(y,2);
    eig_Sw = zeros(2,n+1);
    lambda = eig(Sw0);
    eig_Sw(:,1) = [min(lambda); max(lambda)];
    for j=1:n
        gmp_up.updatePos(x(j), y(:,j));

        lambda = eig(gmp_up.getSigmaW());
        eig_Sw(:,j+1) = [min(lambda); max(lambda)];
    end
    
    W = gmp.W;
    Sw = gmp_up.getSigmaW();

end

function [W, Sw, eig_Sw] = recursive_UD_LS(y, y_old, H, W0, Sw0, R, R_old)

    error('TODO');
    
    if (nargin < 6), R=1; end
    if (nargin < 7), R_old=1; end
    
    W = W0;
    Sw = Sw0;%eye(N_kernels,N_kernels);%Sw;
    n = length(y);
    eig_Sw = zeros(2,n+1);
    eig_Sw(:,1) = [Sw(1,1); Sw(1,1)];
    for i=1:n

        % downdate
        h = H(i,:);
        a = 1/(-R_old + h*Sw*h');
        b = Sw*h'*a;
        W = W + b*(y_old(i) - h*W);
        Sw = Sw - b*h*Sw;

        % update
        h = H(i,:);
        a = 1/(R + h*Sw*h');
        c = h*Sw;
        b = c'*a;
        W = W + b*(y(i) - h*W);
        Sw = Sw - b*c;

        lambda = eig(Sw);
        eig_Sw(:,i+1) = [min(lambda); max(lambda)];
    end

end

%% Plot functions

function fig = plotEigSw(eig_Sw, p)
    
    if (nargin < 2), p = [0.95 0.05]; end
    
    n = size(eig_Sw,2);
    
    fig = figure;
    k = floor(p*n);
    ylabels = {'$\lambda_{min}$', '$\lambda_{max}$'};
    for i=1:2
        ax=subplot(2,1,i); hold on;
        ind = n-k(i):n;
        plot(ind,eig_Sw(i,ind), 'LineWidth',2, 'Color','red', 'LineStyle','-');
        plot(ax.XLim, [0 0], 'LineWidth',1.2, 'Color','magenta', 'LineStyle',':');
        if (i == 2), xlabel('iter \#', 'interpreter','latex', 'fontsize',15); end
        ylabel(ylabels{i}, 'interpreter','latex', 'fontsize',15);
        axis tight;
        set(ax, 'YScale', 'log');
        
        ax.YTick = sort([ax.YTick eig_Sw(i,end)]);
%         ax.YTickLabel = num2str(ax.YTick * double(10^(-ax.YAxis.Exponent)) );
        
    end

end

function fig = plotSw(Sw)
    
    fig = figure;
    imshow(Sw, [min(Sw(:)) max(Sw(:))], 'InitialMagnification', 2000);

end
