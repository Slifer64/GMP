function [Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp = gmp0.deepCopy();
    
    n_dof = length(y0);
    
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);

    gmp_opt = GMP_Opt(gmp);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    n_points = 200;
    % position constr
    gmp_opt.setPosBounds(pos_lim(:,1), pos_lim(:,2), n_points);
    gmp_opt.setPosConstr([],[],[], [0 1], [y0 yg]);
    % velocity constr
    gmp_opt.setVelBounds(vel_lim(:,1), vel_lim(:,2), n_points);
    gmp_opt.setVelConstr([], [], [], [0 1], zeros(n_dof,2));
    % accel constr
    gmp_opt.setAccelBounds(accel_lim(:,1), accel_lim(:,2), n_points);
    gmp_opt.setAccelConstr([], [], [], [0 1], zeros(n_dof,2));
    
    n_points

    % gmp_opt.optimize(100);
    t_start = tic;
%     gmp_opt.setQPsolver(GMP_Opt.MATLAB_QUADPROG);
    gmp_opt.setQPsolver(GMP_Opt.OSQP);
%     gmp_opt.setQPsolver(GMP_Opt.GOLDFARB_IDNANI);
    gmp_opt.optimize2(0:0.01:1);
    elaps_t = toc(t_start);
    fprintf('===> offline-GMP-weights optimization finished! Elaps time: %f ms\n',elaps_t*1000);
    fprintf([gmp_opt.getExitMsg() '\n']);

    [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
    
end
