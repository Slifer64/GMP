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
