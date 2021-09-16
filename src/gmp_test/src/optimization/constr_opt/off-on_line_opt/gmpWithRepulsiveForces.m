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
    dt = 0.001;
    
    y_min = pos_lim(:,1);
    y_max = pos_lim(:,2);
    
    ydot_min = vel_lim(:,1);
    ydot_max = vel_lim(:,2);
    
    yddot_min = accel_lim(:,1);
    yddot_max = accel_lim(:,2);
      
    K = 300;
    D = 80;
    
    ydot_prev = ydot;

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
        
        kp = 0.0001;
        kv = 0.05;
    
        fv = kv * ( 1./(ydot - ydot_max).^3 + 1./(ydot - ydot_min).^3 );
        fp = kp * ( 1./(y - y_max).^3 + 1./(y - y_min).^3 );
        
        z_dot = -K*(y - yd) - D*(z - yd_dot) + yd_ddot + 1*fv + 0*fp;
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

        if (x >= 1.0), break; end

    end
    
    ddP_data = reshape( cell2mat( arrayfun( @(i) [diff(dP_data(i,:))]/dt, 1:n_dof, 'uni',false ) ), size(dP_data,2)-1, size(dP_data,1) )';
    ddP_data = [ddP_data ddP_data(:,end)];
end
