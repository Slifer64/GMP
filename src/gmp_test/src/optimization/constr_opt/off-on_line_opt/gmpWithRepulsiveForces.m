function [Time, P_data, dP_data, ddP_data] = gmpWithRepulsiveForces(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim)
    
    gmp = gmp0.deepCopy();
    
    n_dof = length(y0);
    
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    y = y0;
    ydot = zeros(n_dof,1);
    z = ydot;
    z_dot = zeros(n_dof,1);
    
    gmp.setY0(y0);
    gmp.setGoal(yg);

    t = 0;
    dt = 0.0002;
    
    y_min = pos_lim(:,1);
    y_max = pos_lim(:,2);
    
    ydot_min = vel_lim(:,1);
    ydot_max = vel_lim(:,2);
    
    yddot_min = accel_lim(:,1);
    yddot_max = accel_lim(:,2);
      
    K = 300;
    D = 80;
    
    ydot_prev = ydot;
    
    %% ---------------------------------------------------------
    
    online_plot = false;
    
    if (online_plot)
        fig = figure;
        fig.Position(3:4) = [735 845];
        ax_ = cell(3,1);
        pl_ = cell(3,1);

        plot_step = 10;
        plot_s_i = 0;

        ax_limits = [y_min, ydot_min, yddot_min; y_max, ydot_max, yddot_max];

        [Time0, P0_data, dP0_data, ddP0_data] = getGMPTrajectory(gmp, tau, y0, yg);
        data = {P0_data, dP0_data, ddP0_data};

        y_limits = [y_min, ydot_min, yddot_min; y_max, ydot_max, yddot_max];

        for j=1:3

            ax_{j} = subplot(3,1,j);       
            axes(ax_{j});

            hold on;
            plot([0 tau], ax_limits(1,j)*ones(1,2), 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
            plot([0 tau], ax_limits(2,j)*ones(1,2), 'LineWidth',2, 'LineStyle','--', 'Color','magenta');

            plot(Time0, data{j}, 'LineWidth',2, 'LineStyle',':', 'Color','blue');

            pl_{j} = plot(nan, nan, 'LineWidth',2, 'LineStyle','-', 'Color',[0.85 0.33 0.1]);

            xlim([0 tau]);

        end
    
    end
    
    tic
    
    %% ---------------------------------------------
    while (true)

        x = t/tau;
        x_dot = 1/tau;
        x_ddot = 0;
        
        if (x >= 1)
            x_dot = 0;
        end
        
        yd = gmp.getYd(x);
        yd_dot = gmp.getYdDot(x, x_dot);
        yd_ddot = gmp.getYdDDot(x, x_dot, x_ddot);
        
%         fv = -kv ./ ( (1 - ydot./ydot_max) .* (ydot./ydot_min - 1) );
%         fp = kp ./ ( (1 - y./y_max) .* (y./y_min - 1) );
        
        kp = 0.1; % 0.00001;
        kv = 1; %0.0001;
    
        fv = rep_force(ydot,ydot_max,kv) + rep_force(ydot,ydot_min,kv);
        fp = rep_force(y,y_max,kp) + rep_force(y,y_min,kp);
        
        z_dot = -K*(y - yd) - D*(z - yd_dot) + yd_ddot + 1*fv + 1*fp;
        ydot = z + 0*fp;

        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data ydot];
        
        yddot = (ydot - ydot_prev) / dt;
%         if ( ~isempty(find(yddot < yddot_min)) || ~isempty(find(yddot > yddot_max)) )
%             warning('Acceleration limits violated!');
%             break;
%         end
        
        ydot_prev = ydot;

        t = t + dt;
        y = y + ydot*dt;
        z = z + z_dot*dt;
        
        if (online_plot)
            dat = [y ydot yddot];
            for j=1:3
                pl_{j}.XData = [pl_{j}.XData t];
                pl_{j}.YData = [pl_{j}.YData dat(j)];
            end
            plot_s_i = plot_s_i + 1;
            if (plot_s_i == plot_step)
                plot_s_i = 0;
                drawnow
                pause
            end
        end

        if (x >= 1.0), break; end

    end
    
    fprintf('===> GMP with rep-forces finished! Elaps time: %f ms\n',toc()*1000);
    
    ddP_data = reshape( cell2mat( arrayfun( @(i) [diff(dP_data(i,:))]/dt, 1:n_dof, 'uni',false ) ), size(dP_data,2)-1, size(dP_data,1) )';
    ddP_data = [ddP_data ddP_data(:,end)];
    
end

function f = rep_force(y, y_lim, gain)

    
    d0 = 0.01;
    x = y-y_lim;
    x_norm = abs(y-y_lim);
    
%     scale = 10;
%     a = 20;
%     k = gain * 1 / ( 1 + exp(a*scale*(x_norm-d0)) );
%     f = k * 1./(y - y_lim).^3;
 
    f = 0;
    if (x_norm < d0)
        e = (d0 - x_norm)*x/x_norm;
        psi = (x_norm - d0)^2 / d0^2;
        f = ( 2*gain / ( d0^2*(1-psi) ) ) * log(1/(1-psi))*e;
    end


end
