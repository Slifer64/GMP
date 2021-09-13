function [Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp2 = gmp.deepCopy();
    
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    gmp_opt = GMP_Opt(gmp2);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    n_points = 200;
    % position constr
    gmp_opt.setPosBounds(pos_lim(:,1), pos_lim(:,2), n_points);
    gmp_opt.setPosConstr([],[],[], [0 1], [y0 yg]);
    % velocity constr
    gmp_opt.setVelBounds(vel_lim(:,1), vel_lim(:,2), n_points);
    gmp_opt.setVelConstr([], [], [], [0 1], zeros(3,2));
    % accel constr
    gmp_opt.setAccelBounds(accel_lim(:,1), accel_lim(:,2), n_points);
    gmp_opt.setAccelConstr([], [], [], [0 1], zeros(3,2));
    
    n_points

    % gmp_opt.optimize(100);
    tic
    gmp_opt.optimize2(0:0.01:1);
    elaps_t = toc;
    fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
    fprintf([gmp_opt.getExitMsg() '\n']);
    
    [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
    
end
