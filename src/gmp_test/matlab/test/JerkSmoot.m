classdef JerkSmoot < handle
    
    methods (Access = public, Static)
        
        function run()
            
            model = JerkSmoot();
            
            Time = [];

            x_data = [];
            x_dot_data = [];
            x_ddot_data = [];

            P_data = [];
            P_dot_data = [];
            P_ddot_data = [];
            
            Pd_data = [];
            Pd_dot_data = [];
            Pd_ddot_data = [];

            dt = 0.002;
            t = 0;
            x = 0;
            D = 50;
            K = 250;
            
            model.D = D;
            model.K = K;
            
            dx_ref0 = 1/12;

            P0 = [-0.3023;  -0.0798;   0.6317];  
            Pg = [-0.6040;   0.0251;   0.3000];
            Q0 = [0.0040;   0.8888;   0.4060;  -0.2127];  Q0=Q0/norm(Q0);
            Qg = [0.0097;   0.6711;   0.7411;  -0.0142];  Qg=Qg/norm(Qg);
            
            model.setTargetPose(Pg, Qg);
            model.setInitialPose(P0, Q0);
            
            P = P0;
            P_dot = zeros(3,1);
            P_ddot = zeros(3,1);
            
            Pg_dot = zeros(3,1);
            
            Pg2 = Pg;
            
            x = 0;
            x_dot = 0;
            x_ddot = 0;
             
            model.x = x;
            model.x_dot = x_dot;
            model.x_ddot = x_ddot;
            model.P = P;
            model.P_dot = P_dot;
            model.P_ddot = P_ddot;
            
%             x = 0:0.001:0.1;
%             y = 1 ./ (1 + exp(250*(x-0.05)) );
%             figure;
%             plot(x,y, 'LineWidth',2);
%             return;
            
            count = 1;
            iter = 1;

%             [f, df] = model.goalChangeObjFun(Pg)
%             [c, ceq, dc, dceq] = model.goalChangeObjFunCon(Pg, 0.2)
%             
%             return

            Pg_opt_data = [Pg];
            Pg_data = [Pg];
            obj_fun_data = [];
            exit_flag_data = [];
            
            while (true)
                
                model.setTargetPose(Pg, Qg);
               
                if (x>=1), break; end
                
                if (x >= 0.2 && count==1)
                    count = 2;
                    Pg2 = [-0.664; 0.0205; 0.316]; % Pg + [0.1; -0.15; 0.2]; %
                end
                
                % get phase var
                Pd = model.getRefPos(x);
                x_stop = 1 / ( 1 + exp(250*(norm(P-Pd) - 0.05)) );
                dx_ref = dx_ref0 * x_stop;
                x_ddot = -50*(x_dot - dx_ref);
                
                % update current state
                model.x = x;
                model.x_dot = x_dot;
                model.x_ddot = x_ddot;
                model.P = P;
                model.P_dot = P_dot;
                model.P_ddot = P_ddot;

%                 Pg_dot = 2*(Pg2 - Pg);
                
                if (norm(Pg - Pg2) > 5e-3)
                    
                    c_a = 0.4;
                    
                    Pg_new = model.getNewGoal(model, Pg2, Pg, P0, P, P_dot, P_ddot, x, x_dot, x_ddot, D, K, c_a);
                    ex_flag = 2;
                    fval = 0.5*norm(Pg_new - Pg2)^2;
                    
                    pause
                    
%                     options = optimoptions('fmincon','MaxIterations',100, 'SpecifyObjectiveGradient',true, ...
%                         'SpecifyConstraintGradient',false, 'Display','off');
%                     obj_fun = @(s)goalChangeObjFun(model,s, Pg2);
%                     obj_fun_con = @(s)goalChangeObjFunCon(model,s, Pg2, c_a);
%                     [Pg_new, fval, ex_flag] = fmincon(obj_fun,Pg,[],[],[],[],[],[],obj_fun_con,options);
                    
                    obj_fun_data = [obj_fun_data fval];
                    exit_flag_data = [exit_flag_data ex_flag];
                    
                    Pg_dot = (Pg_new - Pg)/dt;
                    %Pg_dot = 100*(Pg_new - Pg); %/dt;
                    
%                     '-------------------'
%                     iter
%                     fval
                    iter = iter + 1;
                    
                    Pg_data = [Pg_data Pg];
                    Pg_opt_data = [Pg_opt_data Pg_new];
                else
                    %Pg_dot = (Pg2 - Pg)/dt;
                    Pg_dot = 2*(Pg2 - Pg);
                end

                % get reference trajectory
                Pd = model.getRefPos(x);
                Pd_dot = model.getRefVel(x, x_dot);
                Pd_ddot = model.getRefAccel(x, x_dot, x_ddot);
                
                P2_ddot = Pd_ddot - D*(P_dot - Pd_dot) - K*(P - Pd);
                
%                 jerk = P2_ddot - P_ddot;
%                 ac = zeros(3,1);
%                 for i=1:3
%                   da = jerk(i);
%                   if (da < 1), ac(i) = 40;
%                   elseif (da > 3), ac(i) = 40*exp(-3*(3-1));
%                   else, ac(i) = 40*exp(-3*(da-1));
%                   end
%                 end
%                 P_3dot = ac .* jerk;
                P_3dot = (P2_ddot - P_ddot) / dt;
                
                % logging
                Time = [Time t];
                x_data = [x_data x];
                x_dot_data = [x_dot_data x_dot];
                x_ddot_data = [x_ddot_data x_ddot];
                P_data = [P_data P];
                P_dot_data = [P_dot_data P_dot];
                P_ddot_data = [P_ddot_data P_ddot];
                Pd_data = [Pd_data Pd];
                Pd_dot_data = [Pd_dot_data Pd_dot];
                Pd_ddot_data = [Pd_ddot_data Pd_ddot];
                
                % integration
                t = t + dt;
                x = x + x_dot*dt;
                x_dot = x_dot + x_ddot*dt;
                P = P + P_dot*dt;
                P_dot = P_dot + P_ddot*dt;
                P_ddot = P_ddot + P_3dot*dt;
                
                Pg = Pg + Pg_dot*dt;
                
            end
            
            n = length(Time)
            
            
            jerk_data = zeros(3,n-1);
            for i=1:3, jerk_data(i,:) = diff(P_ddot_data(i,:)); end
            max_jerk = 0;
            t_max_jerk = 0;
            max_accel = 0;
            t_max_accel = 0;
            for j=10:n-1
                jerk = norm(jerk_data(:,j));
                accel = norm(P_ddot_data(:,j));
                if (jerk > max_jerk)
                    max_jerk=jerk;
                    t_max_jerk = Time(j);
                end
                if (accel > max_accel)
                    max_accel=accel;
                    t_max_accel = Time(j);
                end
            end
            t_max_jerk
            max_jerk
            t_max_accel
            max_accel
            
            figure;
            Time2 = (0:size(Pg_opt_data,2)-1)*dt;
            for i=1:3
               subplot(3,1,i); hold on;
               plot(Time2, Pg_opt_data(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
               plot(Time2, Pg_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
               scatter(Time2(end), Pg2(i), 'Marker','o', 'LineWidth',2, 'SizeData',100, 'MarkerEdgeColor','red');
            end
            
            figure; hold on;
            plot(obj_fun_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
            for i=1:length(exit_flag_data)
                if (exit_flag_data(i) == 0)
                    scatter(i, obj_fun_data(i), 'Marker','o', 'LineWidth',2, 'SizeData',100, 'MarkerEdgeColor',[0.85 0.33 0.1]);
                elseif (exit_flag_data(i) == -2)
                    scatter(i, obj_fun_data(i), 'Marker','o', 'LineWidth',2, 'SizeData',100, 'MarkerEdgeColor','red');
                end
            end
            
            %% Plot results
            i0 = 4;
            for i=1:3
                figure;
                subplot(3,1,1);
                hold on;
                plot(Time(i0:end), P_data(i,i0:end), 'LineWidth',2.0 , 'Color','blue');
                plot(Time(i0:end), Pd_data(i,i0:end), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
                ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
                legend({'actual','ref'}, 'interpreter','latex', 'fontsize',15);

                axis tight;
                hold off;

                subplot(3,1,2);
                hold on;
                plot(Time(i0:end), P_dot_data(i,i0:end), 'LineWidth',2.0, 'Color','blue');
                plot(Time(i0:end), Pd_dot_data(i,i0:end), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
                ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
                axis tight;
                hold off;

                subplot(3,1,3);
                hold on;
                plot(Time(i0:end), P_ddot_data(i,i0:end), 'LineWidth',2.0, 'Color','blue');
                plot(Time(i0:end), Pd_ddot_data(i,i0:end), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
                ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
                axis tight;
                hold off;
            end
            
            
            
        end
        
       
        function Pg_new = getNewGoal(model, Pg2, Pg, P0, P, P_dot, P_ddot, x, x_dot, x_ddot, D, K, c_a)
            
            % ------------------------------------
            figure;
            for i=1:3
                
                g_vals = -1:0.05:1;
                
                n = length(g_vals);
                
                f_obj = zeros(1,n);
                f_con = zeros(1,n);
                
                gd = Pg2(i);
                
                for j=1:n
                    g = g_vals(j);
                    Pg_ = Pg;
                    Pg_(i) = g;
                    model.setTargetPose(Pg_, model.Qg);

                    Pd = model.getRefPos(x);
                    Pd_dot = model.getRefVel(x, x_dot);
                    Pd_ddot = model.getRefAccel(x, x_dot, x_ddot);

                    P2_ddot = Pd_ddot + D*(Pd_dot - P_dot) + K*(Pd - P);
                    
                    f_obj(j) = abs(g - gd);
                    f_con(j) = abs(P_ddot(i) - P2_ddot(i)) - c_a;
                end
                
                subplot(3,1,i); hold on;
                plot(g_vals, f_obj, 'LineWidth',2, 'Color','blue');
                plot(g_vals, f_con, 'LineWidth',2, 'Color','magenta');
%                 ylim([0 2]);
                
            end
            
            model.setTargetPose(Pg2, model.Qg);
            Pd = model.getRefPos(x);
            Pd_dot = model.getRefVel(x, x_dot);
            Pd_ddot = model.getRefAccel(x, x_dot, x_ddot);
            P2_ddot = Pd_ddot + D*(Pd_dot - P_dot) + K*(Pd - P);
            f_con = abs(P_ddot - P2_ddot) - c_a
            Pg2

            model.setTargetPose(Pg, model.Qg);
            
            % ------------------------------------
            
            Pg_new = zeros(3,1);
            
            model.setTargetPose(Pg2, model.Qg);
            Pd = model.getRefPos(x);
            Pd_dot = model.getRefVel(x, x_dot);
            Pd_ddot = model.getRefAccel(x, x_dot, x_ddot);
            
            P2_ddot = Pd_ddot + D*(Pd_dot - P_dot) + K*(Pd - P);
            
            c1 = (Pd_ddot + D*Pd_dot + K*Pd) ./ (Pg2 - P0);
            c2 = -D*P_dot - K*P;
            a = c1;
            b = -P0.*c1 + c2 - P_ddot;
            
            model.setTargetPose(Pg, model.Qg);

            
            for i=1:3
                g_new = Pg2(i);
                if (abs(P2_ddot(i) - P_ddot(i)) < c_a)
                    Pg_new(i) = g_new;
                else
                    g1 = (g_new + b(i) + c_a) / (1-a(i));
                    g2 = (g_new + b(i) - c_a) / (1-a(i));
                    
                    if (abs(g1-g_new) < abs(g2-g_new)), Pg_new(i) = g1;
                    else, Pg_new(i) = g2;
                    end
                end
            end
            
            f1 = abs(Pg - Pg2)
            f2 = abs(Pg_new - Pg2)
            pause

        end
        
    end
    
    methods (Access = public)

        function this = JerkSmoot()
            
            set_matlab_utils_path();
            this.gmp = GMP_nDoF.importFromFile('data/pih_ur_gmp_target_model.bin', 'pos_');
            
        end

        function p = getRefPos(this, x)

            p = this.t2w_pos(this.gmp.getYd(x));

        end

        function v = getRefVel(this, x, x_dot)

            v = this.t2w_vel(this.gmp.getYdDot(x, x_dot));

        end

        function a = getRefAccel(this, x, x_dot, x_ddot)

            a = this.t2w_vel(this.gmp.getYdDDot(x, x_dot, x_ddot));

        end
        
        function setTargetPose(this, Pg, Qg)

            this.Pg = Pg;
            this.Qg = Qg;

            this.Qt = quatinv(Qg')';
            this.Rt = quat2rotm(this.Qt');
            this.p_o = Pg;

            this.gmp.setGoal(this.w2t_pos(Pg));

        end
        
        
        function setInitialPose(this, P0, Q0)

            this.P0 = P0;
            this.gmp.setY0(this.w2t_pos(P0));

        end

        function p2 = t2w_pos(this, p)

            p2 = this.Rt'*p + this.p_o;

        end

        function v2 = t2w_vel(this, v)

            v2 = this.Rt'*v;
            
        end
        
        function p2 = w2t_pos(this, p)

            p2 = this.Rt*(p - this.p_o);

        end

        function v2 = w2t_vel(this, v)

            v2 = this.Rt*v;
            
        end
        
        function [f, df] = goalChangeObjFun(this, s, Pg2)
            
            Pg = s;

            f = 0.5*sum((Pg - Pg2).^2);
            df = (Pg - Pg2);

        end
        
        
        function [c, ceq, dc, dceq] = goalChangeObjFunCon(this, s, Pg2, c_a)
            
            Pg = s;
            Pg0 = this.Pg;
            
            this.setTargetPose(Pg, this.Qg);
            Pd = this.getRefPos(this.x);
            Pd_dot = this.getRefVel(this.x, this.x_dot);
            Pd_ddot = this.getRefAccel(this.x, this.x_dot, this.x_ddot);
            P2_ddot = Pd_ddot - this.D*(this.P_dot - Pd_dot) - this.K*(this.P - Pd);
            this.setTargetPose(Pg0, this.Qg);

            c = 0.5*( sum( (P2_ddot - this.P_ddot).^2 ) - c_a^2 );
            
            P0 = this.P0;
            Yref0 = ( Pd_ddot + this.D*Pd_dot + this.K*Pd ) ./ (Pg - P0);
            dc = diag(Yref0) * (P2_ddot - this.P_ddot);
            
            ceq = [];
            dceq = [];

        end

    end
    
    
    properties (Access = public)
        
        gmp
        
        Qt
        Rt
        p_o
        
        x
        x_dot
        x_ddot
        
        K
        D
        
        P
        P_dot
        P_ddot
        
        Pg
        P0
        Qg

    end
    
end



