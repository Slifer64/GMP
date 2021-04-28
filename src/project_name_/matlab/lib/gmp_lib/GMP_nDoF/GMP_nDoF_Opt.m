%% N-DoF GMP Optimization class
%

classdef GMP_nDoF_Opt < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_nDoF_Opt(gmp)
                
            this.ex_flag_map = containers.Map('KeyType','double','ValueType','char');
            this.ex_flag_map(-10) = 'Empty... Call ''constrOpt'' first.';
            this.ex_flag_map(1) = 'Function converged to the solution x.';
            this.ex_flag_map(0) = 'Number of iterations exceeded options.MaxIterations.';

            this.ex_flag_map(-2) = 'Problem is infeasible. Or, for interior-point-convex, the step size was smaller than options.StepTolerance, but constraints were not satisfied.';
            this.ex_flag_map(-3) = 'Problem is unbounded';
            this.ex_flag_map(2) = 'Step size was smaller than options.StepTolerance, constraints were satisfied.';
            this.ex_flag_map(-6) = 'Nonconvex problem detected.';
            this.ex_flag_map(-8) = 'Unable to compute a step direction.';
            
            this.exit_msg = '';
                
            this.gmp = gmp;
            
            this.setOptions(1,0,0, 1,1,1);
            
        end
        
        function setOptions(this, opt_pos, opt_vel, opt_accel, pos_obj_w, vel_obj_w, accel_obj_w)
            
            this.opt_pos = opt_pos;
            this.opt_vel = opt_vel;
            this.opt_accel = opt_accel;
            
            this.w_p = pos_obj_w;
            this.w_v = vel_obj_w;
            this.w_a = accel_obj_w;

        end       
        
        %% Trajectory optimization.
        %  @param[in] num_points: Number of discreet points to use in the objective function.
        %  @param[in] tau: Duration of motion.
        %  @param[in] pos_constr: Vector of @GMPConstr position constraints. For no constraints pass '[]'.
        %  @param[in] vel_constr: Vector of @GMPConstr velocity constraints. For no constraints pass '[]'.
        %  @param[in] accel_constr: Vector of @GMPConstr acceleration constraints. For no constraints pass '[]'.
        %  @param[out] success: true if minimum found, false otherwise.
        function success = constrOpt(this, num_points, tau, pos_constr, vel_constr, accel_constr)

            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            success = true;
            this.exit_msg = '';
            
            x_data = linspace(0,1, num_points);
            
            % calculate cost function: J = 0.5w'Hw + f'w
            N = length(x_data);
            
            H = 1e-6*ones(n_ker, n_ker); % for numerical stability
            f = zeros(n_ker, n_dof);
            
            sc = this.gmp.getScaling();
            y0 = this.gmp.getY0();
            y0d = this.gmp.getY0d();
            
            if (this.opt_pos)
                H1 = zeros(n_ker, n_ker);
                f1 = zeros(n_ker, n_dof);
                for i=1:N
                    phi = this.gmp.regressVec(x_data(i));
                    H1 = H1 + phi*phi';
                    yd = (this.gmp.getYd(x_data(i)) - y0) + sc*y0d;
                    f1 = f1 - phi*yd';
                end
                H = H + this.w_p*H1;
                f = f + this.w_p*f1;
            end
  
            if (this.opt_vel)
                H2 = zeros(n_ker, n_ker);
                f2 = zeros(n_ker, 1);
                for i=1:N
                    phi = this.gmp.regressVecDot(x_data(i), 1/tau);
                    H2 = H2 + phi*phi';
                    % TODO: fix wrong scalings...
                    f2 = f2 - phi*this.gmp.getYdDot(x_data(i), 1/tau)';
                end
                H = H + this.w_v*H2;
                f = f + this.w_v*f2;
            end
            
            if (this.opt_accel)
                H3 = zeros(n_ker, n_ker);
                f3 = zeros(n_ker, 1);
                for i=1:N
                    phi = this.gmp.regressVecDDot(x_data(i), 1/tau, 0);
                    H3 = H3 + phi*phi';
                    % TODO: fix wrong scalings...
                    f3 = f3 - phi*this.gmp.getYdDDot(x_data(i), 1/tau, 0)';
                end
                H = H + this.w_a*H3;
                f = f + this.w_a*f3;
            end

            % calculate constraint matrices in canonical form
            [A,b, Aeq,beq] = this.getConstrMat(tau, pos_constr, vel_constr, accel_constr);
            
            opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'Display','off');
            
            % solve optimization problem
            % tic
            inv_ks = this.gmp.getInvScaling();
            W = zeros(n_dof, n_ker);
            for i=1:n_dof
                if (isempty(b)), b_i=[];
                else, b_i = b(:,i);
                end
                
                if (isempty(beq)), beq_i=[];
                else, beq_i = beq(:,i);
                end
                
                [W(i,:), ~, ex_flag] = quadprog(H,f(:,i), A,b_i, Aeq,beq_i, [],[], this.gmp.W(i,:), opt);

                if (ex_flag == 1 || ex_flag == 2)
                    %this.gmp.W(i,:) = w'; %w'/ks(i);
                else
                    success = false;
                    % break; ?
                end
                
                if (~isempty(this.exit_msg)), this.exit_msg = [this.exit_msg '\n']; end
                this.exit_msg = [this.exit_msg 'DoF-' num2str(i) ': ' this.ex_flag_map(ex_flag)];

            end
            
            if (success)
                this.gmp.W = inv_ks*W;
            end
            
            % toc
            
        end
        
        function msg = getExitMsg(this)
           
            msg = this.exit_msg;
            
        end
        
    end
    
    methods (Access = public, Static)
   
        function constr = upperBoundConstr(ti, bound)

            constr = repmat(GMPConstr(), length(ti), 1);
            for i=1:length(ti)
                constr(i) = GMPConstr(ti(i),bound,'<');
            end

        end

        function constr = lowerBoundConstr(ti, bound)

            constr = repmat(GMPConstr(), length(ti), 1);
            for i=1:length(ti)
                constr(i) = GMPConstr(ti(i),bound,'>');
            end

        end

        function constr = thresholdConstr(ti, thres)

            constr = [GMP_nDoF_Opt.upperBoundConstr(ti, thres); GMP_nDoF_Opt.lowerBoundConstr(ti, -thres)];

        end
        
    end
    
    methods (Access = private)

        %% Extracts the constraint matrices in canonical form.
        %  That is: A*w <= b, Aeq*w = beq
        %  Each input argument constraint is of the form struct('t',t_c,'value',val, 'type','=/</>').
        %  @param[in] pos_constr: Vector of position constraints.
        %  @param[in] vel_constr: Vector of velocity constraints.
        %  @param[in] accel_constr: Vector of acceleration constraints.
        %  @param[out] A: Inequality constraint matrix.
        %  @param[out] b: Inequality constraint vector.
        %  @param[out] Aeq: Equality constraint matrix.
        %  @param[out] beq: Equality constraint vector.
        function [A,b, Aeq,beq] = getConstrMat(this, tau, pos_constr, vel_constr, accel_constr)
           
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            
            y0 = this.gmp.getY0();
            y0d = this.gmp.getY0d();
            ks = this.gmp.getScaling();
            
            x_dot = 1/tau;
            
            % position contraints
            if (~isempty(pos_constr))
                for k=1:length(pos_constr)
                    x_c = pos_constr(k).x;
                    val = pos_constr(k).value;
                    type = pos_constr(k).type;
                    
                    phi = this.gmp.regressVec(x_c);                   
                    b_k = ( (val - y0) + ks*y0d )';
                    if (type == '=')
                        Aeq = [Aeq; phi'];
                        beq = [beq; b_k];
                    elseif (type == '<')
                        A = [A; phi'];
                        b = [b; b_k];
                    elseif (type == '>')
                        A = [A; -phi'];
                        b = [b; -b_k];
                    end
                end    
            end
            
            
            % velocity contraints
            if (~isempty(vel_constr))
                for k=1:length(vel_constr)
                    x_c = vel_constr(k).x;
                    val = vel_constr(k).value;
                    type = vel_constr(k).type;
                    
                    phi = this.gmp.regressVecDot(x_c, x_dot);                 
                    b_k = val';
                    
                    if (type == '=')
                        Aeq = [Aeq; phi'];
                        beq = [beq; b_k];
                    elseif (type == '<')
                        A = [A; phi'];
                        b = [b; b_k];
                    elseif (type == '>')
                        A = [A; -phi'];
                        b = [b; -b_k];
                    end
                end  
            end

            % acceleration contraints
            if (~isempty(accel_constr))
                for k=1:length(accel_constr)
                    x_c = accel_constr(k).x;
                    val = accel_constr(k).value;
                    type = accel_constr(k).type;
                    
                    phi = this.gmp.regressVecDDot(x_c, x_dot, 0);                 
                    b_k = val';
                    if (type == '=')
                        Aeq = [Aeq; phi'];
                        beq = [beq; b_k];
                    elseif (type == '<')
                        A = [A; phi'];
                        b = [b; b_k];
                    elseif (type == '>')
                        A = [A; -phi'];
                        b = [b; -b_k];
                    end
                end  
            end
            
        end

    end
    
    properties (Access = private)
        
        ex_flag_map % maps the exit flag value to the msg describing the exit flag
        exit_msg
        
        gmp % n_DoF GMP
        
        opt_pos % flag indicating to optimize position
        opt_vel % flag indicating to optimize velocity
        opt_accel % flag indicating to optimize acceleration
        
        w_p % position objective weight
        w_v % velocity objective weight
        w_a % acceleration objective weight

    end
    
end
