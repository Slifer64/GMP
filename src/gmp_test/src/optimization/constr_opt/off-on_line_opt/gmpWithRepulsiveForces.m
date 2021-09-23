function [Time, P_data, dP_data, ddP_data] = gmpWithRepulsiveForces(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim)
    
    global K D kp kv log_data
    
    K = 300;
    D = 80;
    kp = 0.05;
    kv = 0.5;
    
    gmp = gmp0.deepCopy();
    
    n_dof = length(y0);
    
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    y_min = pos_lim(:,1);
    y_max = pos_lim(:,2);
    
    ydot_min = vel_lim(:,1);
    ydot_max = vel_lim(:,2);
    
    yddot_min = accel_lim(:,1);
    yddot_max = accel_lim(:,2);
    
    tic
    
    %% ---------------------------------------------
        
%     [Time, P_data, dP_data, ddP_data] = simEuler(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max);
    [Time, P_data, dP_data, ddP_data] = simEulerWithOdeStep(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max);
%     [Time, P_data, dP_data, ddP_data] = simOde(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max);
    

    fprintf('===> GMP with rep-forces finished! Elaps time: %f ms\n',toc()*1000);
   
    
end

%% =============================================

function [Time, P_data, dP_data, ddP_data] = simEuler(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max)
        
    global K D kp kv
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];
    
    n_dof = length(y0);
    
    y = y0;
    ydot = zeros(n_dof,1);
    yddot = zeros(n_dof,1);
    ydot_prev = ydot;
    z = ydot;
    z_dot = zeros(n_dof,1);

    t = 0;
    x = 0;
    x_dot = 1/tau;
    xd_dot = 1/tau;
    x_ddot = 0;
    dt = 0.0001;
    
    e_data = [];
    
    while (true)

        x_ddot = -60*(x_dot - xd_dot);

        if (x >= 1)
            x_dot = 0;
            x_ddot = 0;
        end

        yd = gmp.getYd(x);
        yd_dot = gmp.getYdDot(x, x_dot);
        yd_ddot = gmp.getYdDDot(x, x_dot, x_ddot);

%         fv = -kv ./ ( (1 - ydot./ydot_max) .* (ydot./ydot_min - 1) );
%         fp = kp ./ ( (1 - y./y_max) .* (y./y_min - 1) );

        fv = rep_force(ydot,ydot_max,kv) + rep_force(ydot,ydot_min,kv);
        fp = rep_force(y,y_max,kp) + rep_force(y,y_min,kp);

        z_dot = -K*(y - yd) - D*(z - yd_dot) + yd_ddot + 1*fv + 1*fp;
        ydot = z + 0*fp;

        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data ydot];
        ddP_data = [ddP_data yddot];

        yddot = (ydot - ydot_prev) / dt;
%         e = [yddot-yddot_min; yddot_max-yddot];
%         if ( ~isempty(find(e < 0)) )
%             %warning('Acceleration limits violated!');
%             %break;
%             e = max(abs(e(e<0)));  
%             xd_dot = (1/tau) * (1/(1+0.1*e));
%         else
%             xd_dot = 1/tau;
%             e = 0;
%         end
%         e_data = [e_data e];

        ydot_prev = ydot;

        t = t + dt;
        x = x + x_dot*dt;
        x_dot = x_dot + x_ddot*dt;
        y = y + ydot*dt;
        z = z + z_dot*dt;

        if (x >= 1.0), break; end

    end
    
        
end

function [Time, P_data, dP_data, ddP_data] = simOde(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max)

    n_dof = length(y0);
    
    dt = 0.005;
    x0 = 0;
    x0_dot = 1/tau;
    y0 = y0;
    y0_dot = zeros(size(y0));
    
    tspan = 0:dt:tau;
    s0 = [x0; x0_dot; y0; y0_dot];
    opts = odeset('Events',@(t,s)ode_stop_fun(t,s)); % 'OutputFcn',@(t,s,flag)ode_output_fun(t,s,flag)
    [Time,s_data] = ode15s(@(t,s)ode_fun(t,s, gmp, tau, y0, y_min, y_max, ydot_min, ydot_max), tspan, s0, opts);

    s_data = s_data';
    P_data = s_data(3:3+n_dof-1,:);
    dP_data = s_data(3+n_dof:end,:);

    ddP_data = reshape( cell2mat( arrayfun( @(i) [diff(dP_data(i,:))]/dt, 1:n_dof, 'uni',false ) ), size(dP_data,2)-1, size(dP_data,1) )';
    ddP_data = [ddP_data ddP_data(:,end)];

end

function [Time, P_data, dP_data, ddP_data] = simEulerWithOdeStep(gmp, tau, y0, y_min, y_max, ydot_min, ydot_max)
         
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];
    
    n_dof = length(y0);
    
    y = y0;
    ydot = zeros(n_dof,1);
    z = ydot;

    t = 0;
    x = 0;
    x_dot = 1/tau;
    dt = 0.005;
    
    while (true)

        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data ydot];
 
        [~,s_data] = ode15s(@(t,s)ode_fun(t,s, gmp, tau, y0, y_min, y_max, ydot_min, ydot_max), [0 dt], [x; x_dot; y; z]);
        
        t = t + dt;
        [x, x_dot, y, z] = unpackOdeState(s_data(end,:)');
        
        ydot = z;

        if (x >= 1.0), break; end

    end
    
    ddP_data = reshape( cell2mat( arrayfun( @(i) [diff(dP_data(i,:))]/dt, 1:n_dof, 'uni',false ) ), size(dP_data,2)-1, size(dP_data,1) )';
    ddP_data = [ddP_data ddP_data(:,end)];
    
        
end

%% =============================================

function f = rep_force(y, y_lim, gain)

    n_dof = length(y);
    
    d0 = 0.03;
    x = y-y_lim;
    x_norm = abs(y-y_lim);
    
%     scale = 10;
%     a = 20;
%     k = gain * 1 / ( 1 + exp(a*scale*(x_norm-d0)) );
%     f = k * 1./(y - y_lim).^3;
 
    f = zeros(n_dof,1);
    for i=1:n_dof
        if (x_norm(i) < d0)
            e = (d0 - x_norm(i))*x(i)/x_norm(i);
            psi = (x_norm(i) - d0)^2 / d0^2;
            f(i) = ( 2*gain / ( d0^2*(1-psi) ) ) * log(1/(1-psi))*e;
        end
    end

end

%% =============================================

function s_dot = ode_fun(t, s, gmp, tau, y0, y_min, y_max, ydot_min, ydot_max)

    global K D kp kv
    
    xd_dot = 1/tau;
    
    [x, x_dot, y, z] = unpackOdeState(s);
    
    y_dot = z;
    x_dot = x_dot;
    x_ddot = -60*(x_dot - xd_dot);
    
    yd = gmp.getYd(x);
    yd_dot = gmp.getYdDot(x, x_dot);
    yd_ddot = gmp.getYdDDot(x, x_dot, x_ddot);

    fv = rep_force(y_dot,ydot_max,kv) + rep_force(y_dot,ydot_min,kv);
    fp = rep_force(y,y_max,kp) + rep_force(y,y_min,kp);

    z_dot = -K*(y - yd) - D*(z - yd_dot) + yd_ddot + 1*fv + 1*fp;
    y_dot = z + 0*fp;
    
    s_dot = packOdeState(x_dot, x_ddot, y_dot, z_dot);
    
end

function status = ode_output_fun(t,s,flag)

    global log_data
    
    if (strcmpi(flag,'init'))
        log_data.Time = [];
        log_data.P_data = [];
        log_data.dP_data = [];
    end
    
    if (isempty(flag))
        for j=1:length(t)
            [x, x_dot, y, z] = unpackOdeState(s(:,j));
            log_data.Time = [log_data.Time t(j)];
            log_data.P_data = [log_data.P_data y];
            log_data.dP_data = [log_data.dP_data z];
        end
    end
    
    status = 0;

end

function [value,isterminal,direction] = ode_stop_fun(t,s)

    [x, x_dot, y, z] = unpackOdeState(s);
    
    value = 1;
    isterminal = 1;
    direction = 0;

    if (x>=1), value=0; end
    
end

function [x, x_dot, y, z] = unpackOdeState(s)

    x = s(1);
    x_dot = s(2);
    
    n_dof = int32( (length(s)-2)/2 );
    y = s(3:3+n_dof-1);
    z = s(3+n_dof:end);

end

function s = packOdeState(x, x_dot, y, z)

    s = [x; x_dot; y; z];

end

