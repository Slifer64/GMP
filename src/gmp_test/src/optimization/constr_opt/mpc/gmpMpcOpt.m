function [Time, P_data, dP_data, ddP_data] = gmpMpcOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, qp_solver_type)
      
    gmp = gmp0.deepCopy();
   
    n_dof = length(y0);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];
    slack_data = [];
    
    O_ndof = zeros(n_dof,1);

    t = 0;
    dt = 0.005;
    s = t/tau;
    s_dot = 1/tau;
    s_ddot = 0;
    y = y0;
    y_dot = O_ndof;
    y_ddot = O_ndof;
    
    slack_gains = [1e5 100 1];
    
    %% --------  GMP - MPC  --------
    gmp_mpc = GMP_MPC(gmp, 10, 0.1, 30, 1.5, opt_pos, opt_vel, slack_gains); % try 12 , 0.025
    
    gmp_mpc.setPosLimits(pos_lim(:,1), pos_lim(:,2));
    gmp_mpc.setVelLimits(vel_lim(:,1), vel_lim(:,2));
    gmp_mpc.setAccelLimits(accel_lim(:,1), accel_lim(:,2));
    
    gmp_mpc.setPosSlackLimit(5e-3);
    gmp_mpc.setVelSlackLimit(0.05);
    gmp_mpc.setAccelSlackLimit(0.2);
    
%     gmp_mpc.setPosSlackLimit(0);
%     gmp_mpc.setVelSlackLimit(0);
%     gmp_mpc.setAccelSlackLimit(0);
    
    gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
    gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, s_dot, 0, 1e-3*ones(n_dof,1));
    
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    text_prog = ProgressText(40);
    text_prog.init();
    
    t_start = tic;

    %% --------  Simulation loop  --------
    while (true)
        
        %% --------  Stopping criteria  --------
        %if (s > 1.0), break; end
        if (s >= 1 && norm(y_dot)<1e-2 && norm(y-yg)<1e-2 ), break; end

        if (s<=1), text_prog.update(100*t/tau); end

%         if (s >= 1)
%             s = 1;
%             s_dot = 0;
%             s_ddot = 0;
%         end     

        [exit_flag, y, y_dot, y_ddot, slack_var] = gmp_mpc.solve(s, s_dot, s_ddot);
        
        % gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

        if (exit_flag)
            warning(gmp_mpc.getExitMsg());
            text_prog.printInNewLine();
            if (exit_flag < 0), break; end
        end  

        %% --------  Log data  --------
        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data y_dot];
        ddP_data = [ddP_data y_ddot];

        slack_data = [slack_data slack_var];

        %% --------  Numerical integration  --------
        t = t + dt;
        s = s + s_dot*dt;
        s_dot = s_dot + s_ddot*dt;
    
    end

    if (~isempty(slack_data))
        
        pos_slack = slack_gains(1) > 0;
        vel_slack = slack_gains(2) > 0;
        accel_slack = slack_gains(3) > 0;
        
        n_slack = size(slack_data,1);
        y_lb = {};
        if (pos_slack), y_lb = [y_lb, {'pos'}]; end
        if (vel_slack), y_lb = [y_lb, {'vel'}]; end
        if (accel_slack), y_lb = [y_lb, {'accel'}]; end
        figure;
        for i=1:n_slack
            subplot(n_slack,1,i);
            plot(Time, slack_data(i,:), 'LineWidth',2, 'Color','red');
            ylabel(y_lb{i}, 'fontsize',15);
            if (i==1), title('slack variables', 'fontsize',17); end
        end
    end
    
    text_prog.update(100);
    fprintf('\n');
    fprintf('===> GMP-MPC optimization finished! Elaps time: %f ms\n',toc(t_start)*1000);

end
