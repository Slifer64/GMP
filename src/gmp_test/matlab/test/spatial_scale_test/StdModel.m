classdef StdModel < Model
    
    %% ==================  methods  =======================
    
    methods (Access = public)
       
        function this = StdModel()
            
            this = this@Model();
            
        end
        
        
        function [p_err, o_err] = train(this, Time, Pd_data, Qd_data)

            P0 = Pd_data(:,1);
            Pg = Pd_data(:,end);
            Q0 = Qd_data(:,1);
            Qg = Qd_data(:,end);

            this.setInitialPose(P0, Q0);
            this.setTargetPose(Pg, Qg);
            
            [p_err, o_err] = train@Model(this, Time, Pd_data, Qd_data);
                       
        end
        
        
    end
    
    %% ------------- private methods ----------------
    methods (Access = protected)
        
        function setTargetPose(this, Pg, Qg)
            
            this.gmp_p.setGoal(Pg);
            this.gmp_o.setQg(Qg);

        end

        function setInitialPose(this, P0, Q0)

            this.gmp_p.setY0(P0);
            this.gmp_o.setQ0(Q0);

        end

        function Qd = getRefQuat(this, x)

            Qd = this.gmp_o.getQd(x);

        end
        
        function qd = getRefOrient(this, x)

            qd = this.gmp_o.getYd(x); % log(Q*inv(Q0))

        end
        
        function p_w = getRefPos(this, x)

            p_w = this.gmp_p.getYd(x);

        end

        function vel_w = getRefVel(this, x, x_dot)

            vel_w = this.gmp_p.getYdDot(x, x_dot);

        end

        function a_w = calcAccel(this, pos, vel, x, x_dot, x_ddot)

            a_w = this.gmp_p.calcAccel(pos, vel, x, x_dot, x_ddot);

        end

        
    end
    
    
    %% ==================  properties  =======================

    
end