classdef TargetModel < Model
    
    %% ==================  methods  =======================
    
    methods (Access = public)
       
        function this = TargetModel()
            
            this = this@Model();
            
            this.p_o = [0; 0; 0];
            this.Qt = [1; 0; 0; 0];
            this.Rt = quat2rotm(this.Qt');
            
            this.P0_ = [1; 1; 1];
            this.Q0_ = [1; 0; 0; 0];

        end
        
        
        function [p_err, o_err] = train(this, Time, Pd_data, Qd_data)

            P0 = Pd_data(:,1);
            Pg = Pd_data(:,end);
            Q0 = Qd_data(:,1);
            Qg = Qd_data(:,end);

            this.setTargetPose(Pg, Qg);

            n_data = size(Pd_data,2);
            Pd_t_data = zeros(3,n_data);
            for j=1:n_data, Pd_t_data(:,j) = this.w2t_pos(Pd_data(:,j)); end

            this.setInitialPose(P0, Q0);
            this.setTargetPose(Pg, Qg);
            
            [p_err, o_err] = train@Model(this, Time, Pd_t_data, Qd_data);
            
        end
        
        
    end
    
    %% ------------- private methods ----------------
    methods (Access = protected)
        
        function setTargetPose(this, Pg, Qg)

            this.Qt = math_.quatInv(Qg);
            this.Rt = quat2rotm(this.Qt');
            this.p_o = Pg;
            
            this.gmp_p.setGoal(this.w2t_pos(Pg));
            this.setInitialPose(this.P0_, this.Q0_); % reset initial pose based on new target

        end

        function setInitialPose(this, P0, Q0)

            this.P0_ = P0;
            this.Q0_ = Q0;
            this.gmp_p.setY0(this.w2t_pos(P0));

        end

        function p_w = getRefPos(this, x)

            p_w = this.t2w_pos(this.gmp_p.getYd(x));

        end
        
        function Qd = getRefQuat(this, x)

            Qd = this.t2w_quat(this.gmp_o.getQd(x));

        end
        
        function qd = getRefOrient(this, x)

            qd = GMPo.quat2q(this.getRefQuat(x), this.Q0_);

        end

        function vel_w = getRefVel(this, x, x_dot)

            vel_w = this.t2w_vel(this.gmp_p.getYdDot(x, x_dot));

        end

        function a_w = calcAccel(this, pos, vel, x, x_dot, x_ddot)

            a_w = this.t2w_vel(this.gmp_p.calcAccel(this.w2t_pos(pos), this.w2t_vel(vel), x, x_dot, x_ddot));

        end

        function p_t = w2t_pos(this, pos) 

            p_t = this.Rt*(pos - this.p_o);

        end

        function vel_t = w2t_vel(this, vel)

            vel_t = this.Rt*vel;

        end

        function p_w = t2w_pos(this, pos)

            p_w = this.Rt'*pos + this.p_o;

        end

        function vel_w = t2w_vel(this, vel)

            vel_w = this.Rt'*vel;

        end
        
        function Q_t = w2t_quat(this, Q_w)
        
            Q_t = math_.quatProd(this.Qt, Q_w);
          
        end

        function Q_w = t2w_quat(this, Q_t)
        
            Q_w = math_.quatProd(math_.quatInv(this.Qt), Q_t);
            
        end
        
    end
    
    
    %% ==================  properties  =======================
    
    properties (Access = protected)
        
        Qt 
        Rt 
        p_o 
        
        P0_
        Q0_
        
    end
    
    
   
    
end