function [Time, P_data, dP_data, ddP_data] = gmpMpcOptCpp(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, call_mode)
      
    if ( strcmpi(call_mode, 'in') )

        %% Write input data to file
        fid = FileIO('data/gmp_mpc_opt_in.bin', bitor(FileIO.out,FileIO.trunc) );
        fid.write('y0', y0);
        fid.write('yg', yg);
        fid.write('tau', tau);
        fid.write('pos_lim', pos_lim);
        fid.write('vel_lim', vel_lim);
        fid.write('accel_lim', accel_lim);
        gmp_.write(gmp, fid, 'gmp_');
        fid.close();

    elseif ( strcmpi(call_mode, 'out') )

        %% Load output data from file
        fid = FileIO('data/gmp_mpc_opt_out.bin', FileIO.in);
        Time = fid.read('Time');
        P_data = fid.read('P_data');
        dP_data = fid.read('dP_data');
        ddP_data = fid.read('ddP_data');
        slack_data = fid.read('slack_data');
        pos_slack = fid.read('pos_slack');
        vel_slack = fid.read('vel_slack');
        accel_slack = fid.read('accel_slack');
        fid.close();

        if (~isempty(slack_data))
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

    end

end
