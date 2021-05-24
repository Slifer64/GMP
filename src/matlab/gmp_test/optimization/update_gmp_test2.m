clc;
close all;
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

gmp_up = GMP_Update(gmp);
gmp_up.enableSigmawUpdate(true);
gmp_up.initSigmaWfromMsr(Timed/Timed(end));
gmp_up.setMsrNoiseVar(1e-3);

% x = Timed/Timed(end);
% gmp.plotPsi(x);
% gmp_up.plotWeightsCovariance();

% return

%% scale spatio-temporally
kt = 1;
ks = 1;
P0d = gmp.getYd(0); %Pd_data(:,1);
Pgd = gmp.getYd(1); %Pd_data(:,end);
P0 = P0d;
Pg = ks*(Pgd-P0d) + P0;
T = Timed(end)/kt;
Time = Timed/kt;
x_data = Time / T;
x_dot = 1/T;
x_ddot = 0;

gmp.setY0(P0);
gmp.setGoal(Pg);

N = length(x_data);
Pref_data = zeros(n_dof,N);
dPref_data = zeros(n_dof,N);
ddPref_data = zeros(n_dof,N);
P_data = zeros(n_dof,N);
dP_data = zeros(n_dof,N);
ddP_data = zeros(n_dof,N);

Fext_data = zeros(n_dof,N);

dt = Time(2) - Time(1);

P = P0;
dP = zeros(3,1);

% create desired human trajectory
n_data = length(Time);
tm = Time(end)/2;
Ph_data = ks.*(Pd_data - P0d) + P0 + [0.25; 0.36; 0.17]*exp(-0.5*(Time-tm).^2/0.9^2);
dPh_data = zeros(3, n_data);
ddPh_data = zeros(3, n_data);
for i=1:3
    dPh_data(i,:) = [diff(Ph_data(i,:)) 0]/Ts;
    ddPh_data(i,:) = [diff(dPh_data(i,:)) 0]/Ts;
end

% figure;
% for i=1:3
%     subplot(3,1,i); hold on;
%     plot(Time, dPh_data(i,:), 'LineWidth',2);
% end
%
% figure;
% for i=1:3
%     subplot(3,1,i); hold on;
%     plot(Time, Pd_data(i,:), 'LineWidth',2);
%     plot(Time, Ph_data(i,:), 'LineWidth',2, 'LineStyle','--');
% end

up_count = 0;

for j=1:N

    t = Time(j);
    x = x_data(j);

    % get reference
    P_ref = gmp.getYd(x);
    dP_ref = gmp.getYdDot(x, x_dot);
    ddP_ref = gmp.getYdDDot(x, x_dot, x_ddot);

    % update weights
    y = P;
    %y = Ph_data(:,j);
    if (norm(y - P_ref) > 1e-5)% && up_count<1900)

        gmp_up.updatePos(x, y);

        up_count = up_count + 1;
    end

%     Fext = Fext_fun(t);
    Fext = 350*(Ph_data(:,j) - P) + 40*(dPh_data(:,j) - dP);

    ddP = ddP_ref + 50*(dP_ref - dP) + 10*(P_ref - P) + 1*Fext;

    % logging
    Pref_data(:,j) = P_ref;
    dPref_data(:,j) = dP_ref;
    ddPref_data(:,j) = ddP_ref;
    P_data(:,j) = P;
    dP_data(:,j) = dP;
    ddP_data(:,j) = ddP;
    Fext_data(:,j) = Fext;

    % integration
    P = P + dP*dt;
    dP = dP + ddP*dt;

end

figure;
for i=1:3
   subplot(3,1,i);
   plot(Time, Fext_data(i,:), 'LineWidth',2);
   axis tight;
end


P2_data = zeros(size(Pd_data));
dP2_data = zeros(size(dPd_data));
ddP2_data = zeros(size(ddPd_data));
for j=1:length(x_data)
   x = x_data(j);
   P2_data(:,j) = gmp.getYd(x);
   dP2_data(:,j) = gmp.getYdDot(x, x_dot);
   ddP2_data(:,j) = gmp.getYdDDot(x, x_dot, x_ddot);
end
figure;
subplot(3,1,1); hold on;
plot(Time, P2_data(1,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, Pref_data(1,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
plot(Time, Ph_data(1,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
ylim([-1 2]);
subplot(3,1,2); hold on;
plot(Time, dP2_data(1,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, dPref_data(1,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
plot(Time, dPh_data(1,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
ylim([-0.8 0.8]);
subplot(3,1,3); hold on;
plot(Time, ddP2_data(1,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
plot(Time, ddPref_data(1,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
plot(Time, ddPh_data(1,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
ylim([-0.5 0.5]);

return

%% Plot results
fig = figure;
k = [0 3 6];
ax = cell(3,3);
for i=1:3
    k = k + 1;
    ax{i,1} = subplot(3,3,k(1));
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Timed, Pref_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Time, Ph_data(i,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
    plot(Time, Pd_data(i,:), 'LineWidth',2.0 , 'Color',[0 1 0 0.6]);
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    if (i==1), legend({'actual','ref','desired','original'}, 'interpreter','latex', 'fontsize',15); end
    min_p = min(Ph_data(i,:));
    max_p = max(Ph_data(i,:));
    ax{i,1}.YLim = [min_p max_p] + 0.5*[-abs(min_p) abs(max_p)];
    hold off;

    ax{i,2} = subplot(3,3,k(2));
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, dPref_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Time, dPh_data(i,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
    plot(Time, dPd_data(i,:), 'LineWidth',2.0 , 'Color',[0 1 0 0.6]);
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    ax{i,3} = subplot(3,3,k(3));
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, ddPref_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Time, ddPh_data(i,:), 'LineWidth',2.0, 'LineStyle','--', 'Color',[0.85 0.33 0.1 0.6]);
    plot(Time, ddPd_data(i,:), 'LineWidth',2.0 , 'Color',[0 1 0 0.6]);
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end
