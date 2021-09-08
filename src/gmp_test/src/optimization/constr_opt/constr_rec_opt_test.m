clc;
close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
fid = FileIO('data/constr_opt_pos_test_train_data.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

Ts = Timed(2)-Timed(1);


%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

kt = 1.5; % temporal scaling
ks = diag([1 1 1]); % spatial scaling
tau = taud/kt;
y0 = yd0 + 0;
% yg = ks*(ygd - yd0) + y0;
% yg = ygd + [0.1; -0.1; 0.2]; view_ = [171.5301, -2.3630];
yg = ygd + [0.7; -0.7; 0.05];  view_ = [171.9421, -3.0690];


%% ======== Limits ==========
%            lower limit     upper limit
pos_lim = [[-1.2 -1.2 0.2]' [1.2 1.2 0.6]'];
vel_lim = [-0.3*ones(3,1) 0.3*ones(3,1)];  % lower and upper limit, same for all DoFs
accel_lim = [-0.4*ones(3,1) 0.4*ones(3,1)];


% --------- Proportional scaling -----------
gmp.setScaleMethod(TrajScale_Prop(3));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','blue', 'legend','prop', 'plot3D',true, 'plot2D',true);

% % --------- Optimized DMP -> POS -----------
% [Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, true, false);
% data{2} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
%     'color',[0.85 0.33 0.1], 'legend','opt-pos', 'plot3D',true, 'plot2D',true);

% ---------- Online GMP optimization ------------
[Time, P_data, dP_data, ddP_data] = getOnlineOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim);
data{2} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','green', 'legend','onlineOpt-pos', 'plot3D',true, 'plot2D',true);


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
    plot(ax.XLim, [vel_lim(i,1) vel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [vel_lim(i,2) vel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
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
    plot(ax.XLim, [accel_lim(i,1) accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [accel_lim(i,2) accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

end



% ======================================================

function [Time, P_data, dP_data, ddP_data] = getOnlineOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim)
    
    gmp2 = gmp.deepCopy();
    
    %% --------  Init sim  --------
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));

    t = 0;
    dt = 0.002;
    x = t/tau;
    x_dot = 1/tau;
    x_ddot = 0;
    
    %% --------  Init MPC  --------
    N = 8;
    
    K = 300;
    D = 80;
    
    n = 6; % state dim
    m = 3; % control input dim
    
    In = eye(n,n);
    
    % H = eye((n+m)*N, (n+m)*N);
    Qi = blkdiag( 2*eye(3,3) , 1*eye(3,3) );
    Ri = 0*0.01*eye(m,m);

    %H = kron(eye(N,N), blkdiag(Qi,Ri) ); % speye?
    H = blkdiag( kron(speye(N), Qi), kron(speye(N), Ri) );
    
    Ai = (In + [zeros(3,3) eye(3,3); -K*eye(3,3) -D*eye(3,3)]*dt);
    Bi = [zeros(3,3); eye(3,3)]*dt;
    
    Ax = kron(speye(N), speye(n)) + kron(sparse(diag(ones(N-1, 1), -1)), -Ai);
    Bu = kron(speye(N), -Bi);
    Aeq = [Ax, Bu];

    u_min = -inf;
    u_max = inf;
    z_min = [pos_lim(:,1); vel_lim(:,1)];
    z_max = [pos_lim(:,2); vel_lim(:,2)];
    
    X_lb = [repmat(z_min, N,1); repmat(u_min*ones(m,1),N,1)];
    X_ub = [repmat(z_max, N,1); repmat(u_max*ones(m,1),N,1)];
    
    Aineq = speye(N*(n+m));
    
%     solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);
%     X_prev = zeros(N*(n+m),1);
    
    % Create an OSQP object
    prob = osqp;
    
    % - OSQP constraints
    q = zeros(N*(n+m),1);
    beq = zeros( n*N , 1);
    
    A = [Aeq; Aineq];
    lb = [beq; X_lb];
    ub = [beq; X_ub];

    % Setup workspace
    prob.setup(H, q, A, lb, ub, 'warm_start',true, 'verbose',false, 'eps_abs',1e-5, 'eps_rel',1e-5);
    
    use_ud = 0;

    %% ===========  Simulation loop  ===========
    while (true)
        
        %% --------  Stopping criteria  --------
        if (x >= 1.0), break; end
        
        t/tau
        
        if (x >= 1)
            x_dot = 0;
        end

        
        %% --------  calc Beq  --------
        qx = zeros(N*n,1);
        qu = zeros(N*m,1);

        xi = x;
        xi_dot = x_dot;
        xi_ddot = x_ddot;
        
        beq = zeros( n*N , 1);
        z0 = [p; p_dot];
            
        pd_i = gmp2.getYd(xi);
        dpd_i = gmp2.getYdDot(xi, xi_dot);
        ddpd_i = gmp2.getYdDDot(xi, xi_dot, xi_ddot);
            
        for i=0:N-1

            ud_i = K*pd_i + D*dpd_i + ddpd_i;
            
            beq((i*n)+1 : (i+1)*n) = use_ud * Bi*ud_i;
            
            qu((i*m)+1 : (i+1)*m) = zeros(m,1); % -Ri*ud_i;
            
            xi = xi + xi_dot*dt;
            xi_dot = xi_dot + xi_ddot*dt;
            
            pd_i = gmp2.getYd(xi);
            dpd_i = gmp2.getYdDot(xi, xi_dot);
            ddpd_i = gmp2.getYdDDot(xi, xi_dot, xi_ddot);
            
            qx((i*n)+1 : (i+1)*n) = -Qi*[pd_i; dpd_i];

        end
        beq(1:n) = beq(1:n) + Ai*z0; % + Bi*ud_0;
        
        q = [qx; qu];
        
        lb(1:n*N) = beq;
        ub(1:n*N) = beq;
        
        %% --------  calc auxiliary control u  --------
        
        u = 0;

        prob.update('q', q, 'l', lb, 'u', ub);
        
        res = prob.solve();

        if ( res.info.status_val ~= 1)
            res.info
            error(res.info.status);
        end

        u = res.x(N*n+1:N*n+m);
        
%         if (~ isempty( find( Aineq * res.x > X_ub +1e-3 ) ) )
%             ind = find( Aineq * res.x < X_ub );
%             temp = Aineq * res.x - X_lb;
%             temp(temp>0)
%             res.info
%             warning('X_ub Constraints violation!'); 
%         end
%         
%         if (~ isempty( find( Aineq * res.x < X_lb -1e-3 ) ) )
%             temp = Aineq * res.x - X_lb;
%             temp(temp<0)
%             res.info
%             warning('X_lb Constraints violation!'); 
%         end
        
        
%         [X, ~, ex_flag, opt_output] = quadprog(H,q, [],[], Aeq,beq, X_lb,X_ub, X_prev, solver_opt);
% 
%         if (ex_flag == 1 || ex_flag == 2)
%             % success
%         else
%             warning(opt_output.message);
%         end
%         
%         u = X(N*n+1:N*n+m);
%         
%         X_prev = X;

        %% --------  calc u_ref  --------
        p_ref = gmp2.getYd(x);
        dp_ref = gmp2.getYdDot(x, x_dot);
        ddp_ref = gmp2.getYdDDot(x, x_dot, x_ddot);
        
        u_ref = use_ud * ( K*p_ref + D*dp_ref + ddp_ref);
        
        %% --------  Simulate dynamics  --------
        p_ddot = -K*p - D*p_dot + u_ref + u;
        
        z = [p; p_dot];
        z2 = Ai*z + Bi*(u_ref + u);
        
        p_ddot = ( z2(4:6) - z(4:6) ) / dt;
        
        %% --------  Log data  --------
        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];
    
        %% --------  Numerical integration  --------
        t = t + dt;
        x = x + x_dot*dt;
        x_dot = x_dot + x_ddot*dt;
%         p = p + p_dot*dt;
%         p_dot = p_dot + p_ddot*dt;
        p = z2(1:3);
        p_dot = z2(4:6);

    end
    
end


function [Time, P_data, dP_data, ddP_data] = getOnlineOptGMPTrajectory2(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp2 = gmp.deepCopy();
    
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    gmp_opt = GMP_Opt(gmp2);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    % gmp_opt.optimize(100);
%     tic
%     gmp_opt.optimize2(0:0.01:1);
%     elaps_t = toc;
%     fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
%     fprintf([gmp_opt.getExitMsg() '\n']);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));
    
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    t = 0;
    dt = 0.002;
    x = t/tau;
    x_dot = 1/tau;
    x_ddot = 0;
    
    pos_lim = 10 * pos_lim;
    vel_lim = 10 * vel_lim;
    accel_lim = 10 * accel_lim;
    
    n_horiz = 4;
    x_horiz = 0:dt:(n_horiz-1)*dt;
    pos_lb = repmat( pos_lim(:,1), 1,n_horiz);
    pos_ub = repmat( pos_lim(:,2), 1,n_horiz);
    vel_lb = repmat( vel_lim(:,1)*ones(3,1), 1,n_horiz);
    vel_ub = repmat( vel_lim(:,2)*ones(3,1), 1,n_horiz);
    accel_lb = repmat( accel_lim(:,1)*ones(3,1), 1,n_horiz);
    accel_ub = repmat( accel_lim(:,2)*ones(3,1), 1,n_horiz);
    
    x_prev = 0;
    p_prev = gmp2.getYd(x);
    dp_prev = gmp2.getYdDot(x, x_dot);
    ddp_prev = gmp2.getYdDDot(x, x_dot, x_ddot);

    while (true)

        x = t/tau;
        x_dot = 1/tau;
        
        t/tau
        
        if (x >= 1)
            x_dot = 0;
        end
           
        % constraints
        gmp_opt.clearConstr();
        gmp_opt.setPosConstr(x_horiz, pos_lb, pos_ub, x_prev, p_prev);
%         gmp_opt.setVelConstr(x_horiz, vel_lb, vel_ub, x_prev, dp_prev);
%         gmp_opt.setAccelConstr(x_horiz, accel_lb, accel_ub, x_prev, ddp_prev);
%         gmp_opt.setPosConstr(x_horiz, pos_lb, pos_ub, [], []);
%         gmp_opt.setVelConstr(x_horiz, vel_lb, vel_ub, [], []);
%         gmp_opt.setAccelConstr(x_horiz, accel_lb, accel_ub, [], []);
        
        % optimize
        success = gmp_opt.optimize3(x_horiz, gmp);
        if (~success), warning([gmp_opt.getExitMsg() '\n']); end
        
        p_ref = gmp2.getYd(x);
        p_ref_dot = gmp2.getYdDot(x, x_dot);
        p_ref_ddot = gmp2.getYdDDot(x, x_dot, 0);
        
        p = p_ref;
        p_dot = p_ref_dot;
        p_ddot = p_ref_ddot;

        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];

        %p_ddot = p_ref_ddot + 30*(p_ref_dot - p_dot) + 100*(p_ref - p);

        t = t + dt;
%         p = p + p_dot*dt;
%         p_dot = p_dot + p_ddot*dt;
        
        x_horiz = x_horiz + dt/tau;
        
        x_prev = x;
        p_prev = p;
        dp_prev = p_dot;
        ddp_prev = p_ddot;

        if (x >= 1.0), break; end

    end
    
end

function [Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp2 = gmp.deepCopy();
    
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    gmp_opt = GMP_Opt(gmp2);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    n_points = 100;
    % position constr
    gmp_opt.setPosBounds(pos_lim(:,1), pos_lim(:,2), n_points);
    gmp_opt.setPosConstr([],[],[], [0 1], [y0 yg]);
    % velocity constr
    gmp_opt.setVelBounds(vel_lim(:,1), vel_lim(:,2), n_points);
    gmp_opt.setVelConstr([], [], [], [0 1], zeros(3,2));
    % accel constr
    gmp_opt.setAccelBounds(accel_lim(:,1), accel_lim(:,2), n_points);
    gmp_opt.setAccelConstr([], [], [], [0 1], zeros(3,2));

    % gmp_opt.optimize(100);
    tic
    gmp_opt.optimize2(0:0.01:1);
    elaps_t = toc;
    fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
    fprintf([gmp_opt.getExitMsg() '\n']);
    
    [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
    
end

function [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg)

    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));
    
    gmp.setY0(y0);
    gmp.setGoal(yg);

    t = 0;
    dt = 0.002;

    while (true)

        x = t/tau;
        x_dot = 1/tau;
        
        if (x >= 1)
            x_dot = 0;
        end
        
        p_ref = gmp.getYd(x);
        p_ref_dot = gmp.getYdDot(x, x_dot);
        p_ref_ddot = gmp.getYdDDot(x, x_dot, 0);
        
        P = p_ref;
        p_dot = p_ref_dot;
        p_ddot = p_ref_ddot;

        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];

        %p_ddot = p_ref_ddot + 30*(p_ref_dot - p_dot) + 100*(p_ref - p);

        t = t + dt;
        p = p + p_dot*dt;
        p_dot = p_dot + p_ddot*dt;

        if (x >= 1.0), break; end

    end

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
    
    patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z1 z1 z1 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z2 z2 z2 z2], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

%     patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

end
