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
            
            H = 1e-7*eye(n_ker, n_ker); % for numerical stability
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
            
%             zero_tol = 1e-8;
%             
%             A2 = A;
%             A2(abs(A)<zero_tol) = 0;
%             n_A = length(A2(:))
%             nnz_A = length(find(abs(A2)>1e-6))
%             
%             H2 = H;
%             H2(abs(H)<zero_tol) = 0;
%             n_H = length(H2(:))
%             nnz_H = length(find(abs(H2)>1e-6))
%             
%             figure;
%             imshow(H2, 'InitialMagnification', 2000);
%             
%             figure;
%             imshow(A2, 'InitialMagnification', 2000);
%             
%             stop
            
%             fid = FileIO('osqp_problem.bin', bitor(FileIO.out, FileIO.trunc));
%             fid.write('x0', this.gmp.W');
%             fid.write('H',H);
%             fid.write('f', f);
%             fid.write('A_lb',A_lb);
%             fid.write('b_lb', b_lb);
%             fid.write('A_ub',A_ub);
%             fid.write('b_ub', b_ub);
%             fid.write('A_eq', A_eq);
%             fid.write('b_eq', b_eq);
%             fid.write('inv_ks', this.gmp.getInvScaling());
%             fid.close();
%             stop

%             fid = FileIO('qp_solution.bin',FileIO.in);
%             W2 = fid.read('W');
%             
%             A*W2 - b
%             Aeq*W2 - beq
%             
%             pause
            
            % opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'LinearSolver','sparse', 'Display','off');
            opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'Display','off');
            
%             H = sparse(H);
%             A = sparse(A);
%             Aeq = sparse(Aeq);
            
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
        

        function success = constrOpt2(this, num_points, tau, pos_bounds, vel_bounds, accel_bounds)

            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            success = true;
            this.exit_msg = '';
            
            x_data = linspace(0,1, num_points);
            
            % calculate cost function: J = 0.5w'Hw + f'w
            N = length(x_data);
            
            H = 1e-7*eye(n_ker, n_ker); % for numerical stability
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
                    f3 = f3 - phi*this.gmp.getYdDDot(x_data(i), 1/tau, 0)';
                end
                H = H + this.w_a*H3;
                f = f + this.w_a*f3;
            end

            x_dot = 1/tau;
            x_ddot = 0;
            
            y0 = this.gmp.getY0();
            y0d = this.gmp.getY0d();
            ks = this.gmp.getScaling();
            
            % Process inequality constraints
            m = pos_bounds.n_constr + vel_bounds.n_constr + accel_bounds.n_constr;
            A = zeros(m, n_ker);

            for i=1:pos_bounds.n_constr
                pos_bounds.lb(:,i) = ( (pos_bounds.lb(:,i) - y0) + ks*y0d );
                pos_bounds.ub(:,i) = ( (pos_bounds.ub(:,i) - y0) + ks*y0d );
            end
     
            b_lb = [pos_bounds.lb'; vel_bounds.lb'; accel_bounds.lb'];
            b_ub = [pos_bounds.ub'; vel_bounds.ub'; accel_bounds.ub'];
            
            xi = pos_bounds.x;
            k = 0;
            for i=1:length(xi), A(k+i,:) = this.gmp.regressVec(xi(i))'; end
            k = k + length(xi);
            
            xi = vel_bounds.x;
            for i=1:length(xi), A(k+i,:) = this.gmp.regressVecDot(xi(i), x_dot)'; end
            k = k + length(xi);
            
            xi = accel_bounds.x;
            for i=1:length(xi), A(k+i,:) = this.gmp.regressVecDDot(xi(i), x_dot, x_ddot)'; end
            
            A = [-A; A];
            b = [-b_lb; b_ub];
            
            % Process equality constraints
            m_eq = pos_bounds.n_eq_constr + vel_bounds.n_eq_constr + accel_bounds.n_eq_constr;
            A_eq = [];
            b_eq = [];
            if (m_eq ~= 0)
                A_eq = zeros(m_eq, n_ker);

                for i=1:pos_bounds.n_eq_constr
                    pos_bounds.value_eq(:,i) = ( (pos_bounds.value_eq(:,i) - y0) + ks*y0d );
                end

                b_eq = [pos_bounds.value_eq'; vel_bounds.value_eq'; accel_bounds.value_eq'];

                xi = pos_bounds.x_eq;
                k = 0;
                for i=1:length(xi), A_eq(k+i,:) = this.gmp.regressVec(xi(i))'; end
                k = k + length(xi);

                xi = vel_bounds.x_eq;
                for i=1:length(xi), A_eq(k+i,:) = this.gmp.regressVecDot(xi(i), x_dot)'; end
                k = k + length(xi);

                xi = accel_bounds.x_eq;
                for i=1:length(xi), A_eq(k+i,:) = this.gmp.regressVecDDot(xi(i), x_dot, x_ddot)'; end
            end
            
            opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off');
            
            % solve optimization problem
            % tic
            W = zeros(n_dof, n_ker);
            for i=1:n_dof
                
                bi = [];
                beq_i = [];
                if (~isempty(b)), bi = b(:,i); end
                if (~isempty(b_eq)), beq_i = b_eq(:,i); end
                
                [W(i,:), ~, ex_flag] = quadprog(H,f(:,i), A,bi, A_eq,beq_i, [],[], this.gmp.W(i,:), opt);

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
                this.gmp.W = this.gmp.getInvScaling()*W;
            end
            
            % toc
            
        end
        
        
        function setMotionDuration(this, tau)
           
            this.tau = tau;
            
        end
        
        
        function setPosConstr(this, x, lb, ub, xeq, pos_eq)
            
            this.pos_bounds = GMP_nDoF_Opt.boundsStruct(x, lb, ub, xeq, pos_eq);
            
        end
        
        function setVelConstr(this)
            
        end
        
        function setAccelConstr(this)
            
        end
        
        
        function msg = getExitMsg(this)
           
            msg = this.exit_msg;
            
        end

    end
    
    methods (Access = public, Static)
        
        %% Return a struct containing the bounds.
        %  @param[in] x: 1 x N vector of canonical timestamps.
        %  @param[in] lb: dim x N vector with the lower bounds.
        %  @param[in] ub: dim x N vector with the upper bounds.
        %  @return Struct with the bounds.
        function b = boundsStruct(x, lb, ub, x_eq, value_eq)
            
            n_constr = length(x);
            if (size(lb,2) ~= n_constr), error(['Lower bounds must have ' num2str(n_constr) ' columns (constraints).']); end
            if (size(ub,2) ~= n_constr), error(['Upper bounds must have ' num2str(n_constr) ' columns (constraints).']); end
            if (size(lb,1) ~= size(ub,1)), error('Lower and Upper bounds must have the same number of rows'); end
            
            n_eq_constr = length(x_eq);
            if (size(value_eq,2) ~= n_eq_constr), error(['Equality values matrix must have ' num2str(n_eq_constr) ' columns (constraints).']); end
            
            b = struct('n_constr',n_constr, 'x',x, 'lb',lb, 'ub',ub, ...
                'n_eq_constr',n_eq_constr, 'x_eq',x_eq, 'value_eq',value_eq);
            
        end
   
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
