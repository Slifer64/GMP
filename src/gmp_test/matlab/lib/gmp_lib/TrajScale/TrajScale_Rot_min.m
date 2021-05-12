%% TrajScale_Rot_min class

classdef TrajScale_Rot_min < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Rot_min()
    
            this@TrajScale(3);
            
        end

    end
    
    methods (Access = protected)
       
        % ------------------------------------------
        
        function sc = calcScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = axang2rotm([k' theta]);
            end
            
            sc = R*norm(this.Yg - this.Y0)/norm(this.Ygd - this.Y0d);
            
        end
        
        function sc = calcInvScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = axang2rotm([k' theta]);
            end
            sc = R'*norm(this.Ygd - this.Y0d)/norm(this.Yg - this.Y0);
            
        end

        % ------------------------------------------
        
    end
    
end
