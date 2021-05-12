clc;
close all;
clear;

set_matlab_utils_path('../');

% y = [5; 2];
% h = [3 -7 4 6 -3; 1 2 4 5 -1];
%
% [Q,R] = qr(h)
%
% x = h\y


%% Generate data
Ts = 0.002;
Time = 0:Ts:10;
% nominal data
[yd_data, dyd_data, ddyd_data] = get5thOrderPol(0, 0.8, Time);
yd_data = yd_data';
dyd_data = dyd_data';
ddyd_data = ddyd_data';

n_data = length(Time);
tm = Time(end)/2;
% modified data
yh_data = yd_data + 0.25*exp(-0.5*(Time'-tm).^2/0.9^2);
dyh_data = [diff(yh_data')'; 0]/Ts;
ddyh_data = [diff(dyh_data')'; 0]/Ts;

tau = Time(end);
x_data = Time / tau;
x_dot = 1/tau;
x_ddot = 0;


%% Calculate regressors
N_kernels = 20;
gmp = GMP_nDoF(1, N_kernels, 1);
% offline_train_mse = gmp.train('LS', Time/Time(end), Pd_data);
% Wo = gmp.W';

H = zeros(n_data, N_kernels);
H_dot = zeros(n_data, N_kernels);
H_ddot = zeros(n_data, N_kernels);
for i=1:length(x_data)
    x = x_data(i);
    H(i,:) = gmp.regressVec(x)';
    H_dot(i,:) = gmp.regressVecDot(x, x_dot)';
    H_ddot(i,:) = gmp.regressVecDDot(x, x_dot, x_ddot)';
end

%% find the points where Pd and Ph differ
ind2 = findUpdatePoints(yd_data, yh_data, 1e-5);
y2 = yh_data(ind2);
H2 = H(ind2,:);
Time2 = Time(ind2);


%% update Weights with LS
[Wo, Swo] = batch_LS(yd_data, H);
[y0_data, dy0_data, ddy0_data] = getOutput(Wo, H, H_dot, H_ddot);

% [Wo_2, ~] = batch_LS(y0_data, H);
% w_err = norm(Wo-Wo_2)


%% update Weights with RLS
W0 = zeros(N_kernels,1);
Sw0 = 1e6*eye(N_kernels,N_kernels);
R = 1;
[W, Sw, eig_Sw] = recursive_LS(yd_data, H, W0, Sw0, R);

[y_data, dy_data, ddy_data] = getOutput(W, H, H_dot, H_ddot);

plotEigSw(eig_Sw, [0.95 0.05]);

% w_err = norm(W-Wo);
% err_Sw = norm(Sw(:) - Swo(:))
% plotSw(Sw);
% plotSw(Swo);

%% update Weights with recursive U/D LS
% y_old = yd_data(ind2);
y_old = H2*W;
% y2_err = norm(y_old-yd_data(ind2))
[W2, Sw2, eig_Sw2] = recursive_UD_LS(y2, y_old, H2, W, Sw, 1, 1);
% [W2, Sw2, eig_Sw2] = recursive_LS(y2, H2, W, Sw, 1);
% [W2, Sw2, eig_Sw2] = recursive_LS(y2, H2, W, Sw, 1e-3);
[y2_data, dy2_data, ddy2_data] = getOutput(W2, H, H_dot, H_ddot);
plotEigSw(eig_Sw2, [0.95 0.95]);

%% plot trajectories
figure;
subplot(3,1,1); hold on;
plot(Time, y_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, y2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
% plot(Time, yd_data, 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
plot(Time, yh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
% plot(Time2, y2, 'LineWidth',2.0, 'LineStyle',':', 'Color',[0 1 0 0.6]);
ylim([min(yh_data) max(yh_data)]);
subplot(3,1,2); hold on;
plot(Time, dy_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, dy2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
plot(Time, dyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
ylim([min(dyh_data) max(dyh_data)]);
subplot(3,1,3); hold on;
plot(Time, ddy_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, ddy2_data, 'LineWidth',2, 'Color','cyan', 'LineStyle',':');
plot(Time, ddyh_data, 'LineWidth',2.0, 'LineStyle','-', 'Color',[0.85 0.33 0.1 0.6]);
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

function [y_data, dy_data, ddy_data] = getOutput(W, H, H_dot, H_ddot)

    y_data = H*W;
    dy_data = H_dot*W;
    ddy_data = H_ddot*W;

end


%% Update weights functions

function [W, Sw] = batch_LS(y, H)

    W = (H'*H)\H'*y;
    Sw = inv(H'*H);

end

function [W, Sw, eig_Sw] = recursive_LS(y, H, W0, Sw0, R)

    if (nargin < 5), R=1; end
    W = W0;
    Sw = Sw0;
    n = length(y);
    eig_Sw = zeros(2,n+1);
    lambda = eig(Sw);
    eig_Sw(:,1) = [min(lambda); max(lambda)];
    for i=1:n
        h = H(i,:);
        a = 1/(R + h*Sw*h');
        c = h*Sw;
        b = c'*a;
        W = W + b*(y(i) - h*W);
        Sw = Sw - b*c;

        lambda = eig(Sw);
        eig_Sw(:,i+1) = [min(lambda); max(lambda)];
    end
    % W2 = W2 + Sw*H'/(eye(n_data,n_data) + H*Sw*H')*(Y - H*W2);
    % Sw = Sw - P*H'/(eye(n_data,n_data) + H*Sw*H')*H*Sw;

end

function [W, Sw, eig_Sw] = recursive_UD_LS(y, y_old, H, W0, Sw0, R, R_old)

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
