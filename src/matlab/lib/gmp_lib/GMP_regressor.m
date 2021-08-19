%% GMP_regressor classs

classdef GMP_regressor < matlab.mixin.Copyable

    methods (Access = public)
        %% Weighted Sum of Gaussians constructor.
        %  @param[in] N_kernels: The number of kernels.
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std. (optional, default=1.0)
        function this = GMP_regressor(N_kernels, kernel_std_scaling)

            if (nargin < 2), kernel_std_scaling = 1.0; end

            this.setKernels(N_kernels, kernel_std_scaling);
            
        end
        
        %% Set the kernels centers and widths.
        %  @param[in] N_kernels: The number of kernels.
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std.
        function setKernels(this, N_kernels, kernel_std_scaling)
            
            this.c = ((1:N_kernels)-1)'/(N_kernels-1);
            this.h = 1./(kernel_std_scaling*(this.c(2:end)-this.c(1:end-1))).^2;
            this.h = [this.h; this.h(end)];
            
        end

        %% ============================================================
        
        %% Returns the scaled regressor vector ks*phi.
        %  @param[in] x: The phase variable (must be in [0 1]).
        %  @return (scaled) regressor vector.
        function phi = regressVec(this, x)
            
            psi = this.kernelFun(x);
            phi = psi / (sum(psi) + this.zero_tol);

        end
        
        %% Returns the scaled regressor vector 1st time derivative ks*phi_dot.
        %  @param[in] x: The phase variable (must be in [0 1]).
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @return (scaled) regressor vector 1st time derivative.
        function phi_dot = regressVecDot(this, x, x_dot)
            
            psi = this.kernelFun(x);
            psi_dot = this.kernelFunDot(x, x_dot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);

            phi = psi / ( sum(sum_psi) + this.zero_tol );
            phi_dot =  ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this.zero_tol);

        end
        
        %% Returns the scaled regressor vector 2nd time derivative ks*phi_ddot.
        %  @param[in] x: The phase variable (must be in [0 1]).
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @param[in] x_ddot: The phase variable 2nd time derivative.
        %  @return (scaled) regressor vector 2nd time derivative.
        function phi_ddot = regressVecDDot(this, x, x_dot, x_ddot)
            
            psi = this.kernelFun(x);
            psi_dot = this.kernelFunDot(x, x_dot);
            psi_ddot = this.kernelFunDDot(x, x_dot, x_ddot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
            sum_psi_ddot = sum(psi_ddot);

            phi = psi / ( sum(sum_psi) + this.zero_tol );
            phi_dot = ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this.zero_tol);
            phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / ( sum_psi + this.zero_tol);  

        end
        
        %% Returns the scaled regressor vector 3rd time derivative ks*phi_3dot.
        %  @param[in] x: The phase variable (must be in [0 1]).
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @param[in] x_ddot: The phase variable 2nd time derivative.
        %  @param[in] x_3dot: The phase variable 3rd time derivative.
        %  @return (scaled) regressor vector 3rd time derivative.
        function phi_3dot = regressVec3Dot(this, x, dx, ddx, d3x)
            
            psi = this.kernelFun(x);
            psi_dot = this.kernelFunDot(x, dx);
            psi_ddot = this.kernelFunDDot(x, dx, ddx);
            psi_3dot = this.kernelFun3Dot(x, dx, ddx, d3x);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
            sum_psi_ddot = sum(psi_ddot);
            sum_psi_3dot = sum(psi_3dot);

            phi = psi / ( sum(sum_psi) + this.zero_tol );
            phi_dot = ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this.zero_tol);
            phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / ( sum_psi + this.zero_tol);
            phi_3dot = (psi_3dot - 3*phi_ddot*sum_psi_dot - 3*phi_dot*sum_psi_ddot - phi*sum_psi_3dot) / ( sum_psi + this.zero_tol);

        end


        %% =============================================================

        function plotPsi(this, x)
            
            Psi = this.kernelFun(x);
            figure;
            hold on;
            for i=1:length(this.c)
               plot(x, Psi(i,:), 'LineWidth',2); 
            end
            hold off

        end
        
        
    end
    
    methods (Access = protected)
        

        %% =============================================================
        
        %% Returns a column vector with the values of the kernel functions.
        %  @param[in] x: The phase variable.
        %  @return: Column vector with the values of the kernel functions.
        function psi = kernelFun(this, x)

            n = length(x);
            psi = zeros(length(this.c), n);
            
            for j=1:n
                psi(:,j) = exp(-this.h.*((x(j)-this.c).^2));
            end 

        end
        
        function psi_dot = kernelFunDot(this, x, x_dot)

            n = length(x);
            psi = this.kernelFun(x);
            psi_dot = zeros(length(this.c), n);
            
            for j=1:n
                a = (x(j)-this.c)*x_dot(j);
                psi_dot(:,j) = -2*this.h.*( psi(:,j).*a);
            end 

        end
        
        function psi_ddot = kernelFunDDot(this, x, x_dot, x_ddot)

            n = length(x);
            psi = this.kernelFun(x);
            psi_dot = this.kernelFunDot(x,x_dot);
            psi_ddot = zeros(length(this.c), n);

            for j=1:n
                a = (x(j)-this.c)*x_dot(j);
                a_dot = (x(j)-this.c)*x_ddot(j) + x_dot(j)^2;
                psi_ddot(:,j) = -2*this.h.*( psi_dot(:,j).*a + psi(:,j).*a_dot ); 
            end 

        end
        
        function psi_3dot = kernelFun3Dot(this, x, x_dot, x_ddot, x_3dot)

            n = length(x);
            psi = this.kernelFun(x);
            psi_dot = this.kernelFunDot(x,x_dot);
            psi_ddot = this.kernelFunDDot(x,x_dot, x_ddot);
            psi_3dot = zeros(length(this.c), n);
            
            for j=1:n
                a = (x(j)-this.c)*x_dot(j);
                a_dot = (x(j)-this.c)*x_ddot(j) + x_dot(j)^2;
                a_ddot = (x(j)-this.c)*x_3dot(j) + 3*x_dot(j)*x_ddot(j);
                psi_3dot(:,j) = -2*this.h.*( psi_ddot(:,j).*a + 2*psi_dot(:,j).*a_dot + psi(:,j).*a_ddot ); 
            end 

        end
        
        
    end
    
    properties (Access = public)

        c % N_kernels x 1 vector with the kernels' centers
        h % N_kernels x 1 vector with the kernels' inverse width
        
    end
    
    properties (Constant, Access = private)
        
        zero_tol = 1e-200 % small value used to avoid divisions with very small numbers
        
    end
    
end
