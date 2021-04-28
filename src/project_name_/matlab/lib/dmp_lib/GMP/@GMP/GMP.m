%% GMP class
%  Generalized movement primitive.
%

classdef GMP < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] N_kernels: number of kernels
        %  @param[in] D: damping.
        %  @param[in] K: stiffness.
        %  @param[in] kernels_std_scaling: Scaling for std of kernels (optional, default=1).
        function this = GMP(N_kernels, D, K, kernels_std_scaling)
                
            if (nargin < 4), kernels_std_scaling = 1.0; end
            
            this.D = D;
            this.K = K;
            this.wsog = WSoG(N_kernels, kernels_std_scaling);

            this.setY0(0);
            this.setGoal(1);
            this.y_dot = 0;
            this.z_dot = 0;
            
        end

        
        %% Trains the GMP.
        %  @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
        %  @param[in] Time: Row vector with the timestamps of the training data points.
        %  @param[in] yd_data: Row vector with the desired position.
        %  @param[out] train_error: The training error expressed as the mse error.
        function train_error = train(this, train_method, Time, yd_data)

            x = Time / Time(end);
            train_error = this.wsog.train(train_method, x, yd_data);
            
        end

        
        %% Returns the derivatives of the GMP states.
        %  @param[in] s: Vector with the phase variable state, i.e. s = [x; x_dot; x_ddot].
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] z: 'z' state of the GMP.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        function update(this, s, y, z, y_c, z_c)

            if (nargin < 5), y_c=0; end
            if (nargin < 6), z_c=0; end
            
            this.z_dot = ( this.goalAttractor(y, z) + this.shapeAttractor(s) + z_c);
            this.y_dot = z + y_c;

        end

        
        %% Returns the 'y' state time derivative.
        %  Call after @update.
        %  @return: time derivative of the 'y' state.
        function y_dot = getYdot(this)
            y_dot = this.y_dot; 
        end
        
        
        %% Returns the 'z' state time derivative.
        %  Call after @update.
        %  @return: time derivative of 'z' state.
        function z_dot = getZdot(this)
            z_dot = this.z_dot; 
        end
        
        
        %% Returns the GMP's acceleration.
        %  Call after @update.
        %  @param[in] yc_dot: time derivative of 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = getYddot(this, yc_dot)
            
            if (nargin < 2), yc_dot = 0; end
            y_ddot = this.getZdot() + yc_dot;

        end
        
        
        %% Calculates the GMP's acceleration.
        %  @param[in] s: Vector with the phase variable state, i.e. s = [x; x_dot; x_ddot].
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] y_dot: time derivative of 'y' state.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        %  @param[in] yc_dot: time derivative of the 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = calcYddot(this, s, y, y_dot, yc, zc, yc_dot)

            if (nargin < 5), yc = 0; end
            if (nargin < 6), zc = 0; end
            if (nargin < 7), yc_dot = 0; end

            z = y_dot - yc;
            z_dot = ( this.goalAttractor(y, z) + this.shapeAttractor(s) + zc);

            y_ddot = (z_dot + yc_dot);

        end
        
        
        %% Returns the number of kernels.
        %  @return: number of kernels.
        function n_ker = numOfKernels(this)
            
            n_ker = this.wsog.numOfKernels();
            
        end
        
        function s = getSpatialScaling(this), s = this.wsog.getSpatialScaling(); end

        
        %% Sets the initial position.
        %  @param[in] y0: initial position.
        function setY0(this, y0)
            
            this.wsog.setStartValue(y0); 
            
        end
        
        
        %% Set goal position.
        %  @param[in] g: goal position.
        function setGoal(this, g)
            
            this.g = g;
            this.wsog.setFinalValue(g); 
            
        end

        
        %% Returns a deep copy of this object.
        %  @return: deep copy of this object.
        function cp_obj = deepCopy(this)
            
            % Make a shallow copy of all properties
            cp_obj = this.copy();
            % Make a deep copy of the pointers
            cp_obj.wsog = this.wsog.deepCopy();

        end
          

        %% Returns the scaled desired position.
        %  @param[in] x: phase variable.
        %  @return: scaled desired position.
        function p_ref = getYd(this, x)
            
            p_ref = this.wsog.output(x);
            
        end
        
        
        %% Returns the scaled desired velocity.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @return: scaled desired velocity.
        function p_ref_dot = getYdDot(this, x, x_dot)
            
            p_ref_dot = this.wsog.outputDot(x, x_dot);
            
        end
        
        
        %% Returns the scaled desired acceleration.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @param[in] x_ddot: 2nd time derivative of the phase variable.
        %  @return: scaled desired acceleration.
        function p_ref_ddot = getYdDDot(this, x, x_dot, x_ddot)
            
            if (nargin < 4), x_ddot = 0; end
            
            p_ref_ddot = this.wsog.outputDDot(x, x_dot, x_ddot);
            
        end

        
        %% Constrained optimization GMP.
        %  Returns a GMP object optimized to minimize the error from the
        %  unconstrained trajectory, subject to the input constraints. 
        %  @param[in] T: Time duration of the motion.
        %  @param[in] pos_constr: Vector of @GMPConstr position constraints. For no constraints pass '[]'.
        %  @param[in] vel_constr: Vector of @GMPConstr velocity constraints. For no constraints pass '[]'.
        %  @param[in] accel_constr: Vector of @GMPConstr acceleration constraints. For no constraints pass '[]'.
        %  @param[in] opt_set: Object of type @GMPOptSet for setting optimization options.
        %  @param[in] N: number of data points generated for the optimization.
        %  @return: GMP object that produces a trajectory that satisfies the constraints.
        function gmp_opt = constrOpt(this, T, pos_constr, vel_constr, accel_constr, opt_set, N)
        
            ref_data = struct('pos',[], 'vel',[], 'accel',[], 'w_p',[], 'w_v',[], 'w_a',[]);
            
            x = (0:(N-1))/(N-1);
            x_dot = 1/T;
            x_ddot = 0;

            if (opt_set.opt_pos)
                ref_data.pos = zeros(1, N);
                for i=1:N, ref_data.pos(i) = this.getYd(x(i)); end
            end
            
            if (opt_set.opt_vel)
                ref_data.vel = zeros(1, N);
                for i=1:N, ref_data.vel(i) = this.getYdDot(x(i), x_dot); end
            end
            
            if (opt_set.opt_accel)
                ref_data.accel = zeros(1, N);
                for i=1:N, ref_data.accel(i) = this.getYdDDot(x(i), x_dot, x_ddot); end
            end
            
            gmp_opt = this.deepCopy();
            gmp_opt.wsog.constrOpt(x, ref_data, T, pos_constr, vel_constr, accel_constr, opt_set);
            
        end
            
        
        %% Updates the weights so that the generated trajectory passes from the given points.
        %  @param[in] s: Matrix of phase variable states, i.e. s = [x; x_dot; x_ddot].
        %  @param[in] z: Row vector with the desired value for each timestamp.
        %  @param[in] type: Row vector with the type of each point (GMP_UPDATE_TYPE).
        %  @param[in] z_var: Row vector with the variance of each point (optional, default = 1e-3).
        function updateWeights(this, s, z, type, z_var)
            
            if (nargin < 5), z_var = 1e-3; end
            
            n = size(s, 2);
            m = size(s,1);
            
            x = s(1,:);
            x_dot = s(2, :);
            x_ddot = zeros(1,n);
            if (m > 2), x_ddot = s(3, :); end
            
            if (isscalar(z_var)), z_var = z_var*ones(n,1); end
            
            k = this.wsog.numOfKernels();
            
            H = zeros(n, k);
            z_hat = zeros(1,n);
            
            for i=1:n
                if (type(i) == GMP_UPDATE_TYPE.POS)
                    Hi = this.wsog.regressVec(x(i))';
                    z_hat(i) = this.wsog.output(x(i));
                elseif (type(i) == GMP_UPDATE_TYPE.VEL)
                    Hi = this.wsog.regressVecDot(x(i), x_dot(i))';
                    z_hat(i) = this.wsog.outputDot(x(i), x_dot(i));
                elseif (type(i) == GMP_UPDATE_TYPE.ACCEL)
                    Hi = this.wsog.regressVecDDot(x(i), x_dot(i), x_ddot(i))';
                    z_hat(i) = this.wsog.outputDDot(x(i), x_dot(i), x_ddot(i));
                end
                H(i,:) = Hi;
            end
            
            e = (z - z_hat)';
            this.wsog.updateWeights(H, e, diag(z_var));
            
        end
        
        %% Export the GMP model to a file.
        % @param[in] filename: The name of the file.
        function exportToFile(this, filename)
            
            fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
            this.writeToFile(fid, '');
            
        end
        
        %% Write the GMP model to a file.
        % @param[in] fid: Object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function writeToFile(this, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            fid.write([prefix 'D'], this.D);
            fid.write([prefix 'K'], this.K);
            fid.write([prefix 'N_kernels'], this.wsog.N_kernels);
            fid.write([prefix 'w'], this.wsog.w);
            fid.write([prefix 'c'], this.wsog.c);
            fid.write([prefix 'h'], this.wsog.h);
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] fid: Object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function readFromFile(this, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            this.D = fid.read([prefix 'D']);
            this.K = fid.read([prefix 'K']);
            this.wsog.N_kernels = fid.read([prefix 'N_kernels']);
            this.wsog.w = fid.read([prefix 'w']);
            this.wsog.c = fid.read([prefix 'c']);
            this.wsog.h = fid.read([prefix 'h']);

            this.wsog.f0_d = dot(this.wsog.regressVec(0),this.wsog.w);
            this.wsog.fg_d = dot(this.wsog.regressVec(1),this.wsog.w);
            this.wsog.setStartValue(this.wsog.f0_d);
            this.wsog.setFinalValue(this.wsog.fg_d);
            
        end
  
    end
    
    methods (Static, Access = public)
        
        %% Import a GMP model from a file.
        % @param[in] filename: The name of the file.
        function gmp = importFromFile(filename)
            
            gmp = GMP(2, 1, 1, 1);
            fid = FileIO(filename, FileIO.in);
            gmp.readFromFile(fid, '');
  
        end
        
    end
    
    methods (Access = {?GMP_nDoF, ?GMPo} )
                
        
        %% Returns the goal attractor.
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] z: 'z' state of the GMP.
        %  @return: goal attractor.
        function goal_attr = goalAttractor(this, y, z)

            goal_attr = this.K*(this.g-y)- this.D*z;

        end
        
        
        %% Returns the shape attractor.
        %  @param[in] s: Vector with the phase variable state, i.e. s = [x; x_dot; x_ddot].
        %  @return: shape attractor.
        function shape_attr = shapeAttractor(this, s)
            
            x = s(1);
            x_dot = s(2);
            x_ddot = 0;
            if (length(s) > 2), x_ddot = s(3); end
            
            yd = this.getYd(x);
            yd_dot = this.getYdDot(x, x_dot);
            yd_ddot = this.getYdDDot(x, x_dot, x_ddot);
            
            f = yd_ddot + this.D*yd_dot + this.K*(yd - this.g);
            
            sAttrGat = 1;
            shape_attr = sAttrGat * f;
            
        end
        

    end
    
    properties (Access = public)
     
        D % damping
        K % stiffness
        
        wsog % WSoG object for encoding the desired position

    end
    
    
    properties (Access = protected)
        
        %% output state
        y_dot % position (y-state) derivative 
        z_dot % scaled velocity derivative(z-state) derivative
        
        g % target/goal

    end
    
end
