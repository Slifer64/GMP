classdef Model < handle
    
    %% ==================  methods  =======================
    
    methods (Access = public)
       
        function this = Model()
            
            this.gmp_p = GMP_nDoF(3, 30, 40, 150, 1);
            this.gmp_o = GMPo(30, 10, 2, 1);
            
        end
        
        
        function [p_err, o_err] = train(this, Time, Pd_data, Qd_data)
            
            p_err = this.gmp_p.train('LS', Time, Pd_data);
            o_err = this.gmp_o.train('LS', Time, Qd_data);
            
        end
        
        
        function simulate(this, P0, Q0, Pg, Qg, tau, x)
               
            %x_dot = 1/tau;
            n_data = length(x);
            this.P_data = zeros(3, n_data);
            %this.dP_data = zeros(3, n_data);
            this.Q_data = zeros(4, n_data);
            this.q_data = zeros(3, n_data);
 
            this.setInitialPose(P0, Q0);
            this.setTargetPose(Pg, Qg);

            for j=1:n_data
                %x = Time(j)/tau;
                this.P_data(:,j) = this.getRefPos(x(j));
                this.Q_data(:,j) = this.getRefQuat(x(j));
                this.q_data(:,j) = this.getRefOrient(x(j));
                %this.dP_data(:,j) = this.getRefVel(x(j), x_dot);
            end
            
        end
        
    end
    
    %% ------------- private methods ----------------
    methods (Abstract, Access = protected)
        
        setTargetPose(this, Pg, Qg)

        setInitialPose(this, P0, Q0)

        Qd = getRefQuat(this, x)
        
        qd = getRefOrient(this, x)
        
        pd = getRefPos(this, x)

        pd_dot = getRefVel(this, x, x_dot)

        pd_ddot = calcAccel(this, pos, vel, x, x_dot, x_ddot)

    end
    
    
    %% ==================  properties  =======================
    
    
    properties (Access = public)

        P_data
        Q_data
        q_data
        
    end
    
    properties (Access = protected)
        
        gmp_p
        gmp_o
        
    end
    
    
   
    
end