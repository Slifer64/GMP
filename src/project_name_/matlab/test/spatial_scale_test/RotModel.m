classdef RotModel < Model
    
    %% ==================  methods  =======================
    
    methods (Access = public)
       
        function this = RotModel()
            
            this = this@Model();
            
        end
        
        
        function [p_err, o_err] = train(this, Time, Pd_data, Qd_data)

            this.P0d = Pd_data(:,1);
            this.Pgd = Pd_data(:,end);
            this.Q0d = Qd_data(:,1);
            this.Qgd = Qd_data(:,end);

            [p_err, o_err] = train@Model(this, Time, Pd_data, Qd_data);
            
        end
        
    end
    
    %% ------------- private methods ----------------
    methods (Access = protected)
        
        function setTargetPose(this, Pg, Qg)
            
            this.Pg = Pg;
            this.Qg = Qg;

        end

        function setInitialPose(this, P0, Q0)

            this.P0 = P0;
            this.Q0 = Q0;

        end

        function p = getRefPos(this, x)
            
            pd = this.gmp_p.getYd(x);
            p = this.posTransMat()*(pd-this.P0d) + this.P0;

        end

        function Qd = getRefQuat(this, x)

            Qd = GMPo.q2quat(this.getRefOrient(x), this.Q0);
            
        end
        
        function qd = getRefOrient(this, x)

            qd = this.rotTransMat()*this.gmp_o.getYd(x); % log(Q*inv(Q0))

        end
        
            
        function vel = getRefVel(this, x, x_dot)

            pd_dot = this.gmp_p.getYdDot(x, x_dot);
            vel = this.posTransMat() * pd_dot;
            
        end

        function accel = calcAccel(this, pos, vel, x, x_dot, x_ddot)

            pd_ddot = this.gmp_p.calcAccel(pos, vel, x, x_dot, x_ddot);
            accel = this.posTransMat() * pd_ddot;

        end
        
        function Sg = posTransMat(this)
           
            p = this.Pg - this.P0;
            pd = this.Pgd - this.P0d;
            
            s = norm(p) / norm(pd);
            n = p / norm(p);
            nd = pd / norm(pd);
            
            k = cross(nd,n);
            theta = acos(dot(n,nd));
            R = axang2rotm([k' theta]);
            
            Sg = s*R;
            
        end
        
        function Sg = rotTransMat(this)
           
            q = GMPo.quat2q(this.Qg, this.Q0);
            qd = GMPo.quat2q(this.Qgd, this.Q0d);
            
            s = norm(q) / norm(qd);
            n = q / norm(q);
            nd = qd / norm(qd);
            
            k = cross(nd,n);
            theta = acos(dot(n,nd));
            R = axang2rotm([k' theta]);
            
            Sg = s*R;
            
        end

        
    end
    
    
    %% ==================  properties  =======================
    
    properties (Access = protected)
        
        P0
        Q0
        Pg
        Qg
        
        P0d
        Q0d
        Pgd
        Qgd
        
        
    end
    
    
   
    
end