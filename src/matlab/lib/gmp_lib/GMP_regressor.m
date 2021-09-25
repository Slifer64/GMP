%% GMP_regressor class
%  This class implements basic functionalities for a regressor vector
%  constructed as a normalized sum of Gaussian kernel (basis) functions. 
%  One can set the number of the kernels, their width and obtain the 
%  regressor vector and its 1st, 2nd and 3rd derivatives.
%  The kernels are equally spaced in [0 1] and take a one-dimensional
%  input, say 'x', for example:
%  phi = regressVec(x)
%  phi_dot = regressVecDot(x, x_dot)
%  phi_ddot = regressVecDDot(x, x_dot, x_ddot)
%  phi_3dot = regressVec3Dot(x, x_dot, x_ddot, x_3dot)
%

classdef GMP_regressor < matlab.mixin.Copyable
    
    %% ====================================================
    %% ===============  Public functions  =================
    %% ====================================================
    
    methods (Access = public)
        %% Constructor. Initialize with the number of kernels and their width.
        %  @param[in] N_kernels: The number of kernels (must be > 1).
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std. (optional, default=1.0)
        function this = GMP_regressor(N_kernels, kernel_std_scaling)

            if (nargin < 2), kernel_std_scaling = 1.0; end

            this.kernel_fun_ptr = @this.kernelFun;

            this.setKernels(N_kernels, kernel_std_scaling);
            
        end
 
        %% Set the kernels centers and widths.
        %  @param[in] N_kernels: The number of kernels.
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std.
        function setKernels(this, N_kernels, kernel_std_scaling)
            
            if (N_kernels < 2), error('[GMP_regressor::setKernels]: At least 2 kernels are required!'); end
            
            this.c = ((1:N_kernels)-1)'/(N_kernels-1);
            
            this.kernel_std_scaling = kernel_std_scaling;
            this.kernel_std = kernel_std_scaling*(this.c(2:end)-this.c(1:end-1));

            % this.h = this.kernel_std * ones(N_kernels,1);
            this.h = 1./(this.kernel_std).^2;
            this.h = [this.h; this.h(end)];
            
            % find range of x outside which psi = 0 due to finite numerical precision
            r_min = realmin;
            this.x_min = this.c(1) - sqrt( -log(r_min) / this.h(1) );
            this.x_max = this.c(end) + sqrt( -log(r_min) / this.h(end) );
            
        end
        
        %% Enable/disable kernels truncation. 
        %  @param[in] set: enable/disable kernels truncation.
        %  @param[in] zero_tol: threshold below which the activation of a kernel is set to zero (optional, default=1e-8).
        function setTruncatedKernels(this, set, zero_tol)
           
            if (nargin < 3), zero_tol = 1e-8; end
            
            % find range of x outside which psi = 0 due to finite numerical precision
            r_min = zero_tol;
            this.x_min = this.c(1) - sqrt( -log(r_min) / this.h(1) );
            this.x_max = this.c(end) + sqrt( -log(r_min) / this.h(end) );
            
            if (set), this.kernel_fun_ptr = @(x)truncGaussKernel(this, x, zero_tol);
            else, this.kernel_fun_ptr = @this.kernelFun;
            end
        end

        %% ============================================================
        
        %% Returns the regressor vector.
        %  @param[in] x: The phase variable.
        %  @return regressor vector.
        function phi = regressVec(this, x)
            
            % take appropriate actions when x causes phi = 0 due to finite
            % numerical precision.
            if (x < this.x_min)
                psi = zeros(length(this.c),1);
                psi(1) = 1;
            elseif (x > this.x_max)
                psi = zeros(length(this.c),1);
                psi(end) = 1;
            else
                psi = this.kernel_fun_ptr(x);
            end
            
            phi = psi / sum(psi);% + this.zero_tol);

        end
        
        %% Returns the regressor vector 1st time derivative phi_dot.
        %  @param[in] x: The phase variable.
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @return regressor vector 1st time derivative.
        function phi_dot = regressVecDot(this, x, x_dot)
            
            if (x < this.x_min || x > this.x_max)
                phi_dot = zeros(length(this.c),1);
                return;
            end
            
            psi = this.kernel_fun_ptr(x);
            psi_dot = this.kernelFunDot(x, x_dot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);

            phi = psi / sum_psi; % + this.zero_tol );
            phi_dot =  ( psi_dot - phi*sum_psi_dot ) / sum_psi; % + this.zero_tol);

        end
        
        %% Returns the regressor vector 2nd time derivative ks*phi_ddot.
        %  @param[in] x: The phase variable.
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @param[in] x_ddot: The phase variable 2nd time derivative.
        %  @return regressor vector 2nd time derivative.
        function phi_ddot = regressVecDDot(this, x, x_dot, x_ddot)
            
            if (x < this.x_min || x > this.x_max)
                phi_ddot = zeros(length(this.c),1);
                return;
            end
            
            psi = this.kernel_fun_ptr(x);
            psi_dot = this.kernelFunDot(x, x_dot);
            psi_ddot = this.kernelFunDDot(x, x_dot, x_ddot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
            sum_psi_ddot = sum(psi_ddot);

            phi = psi / sum_psi; % + this.zero_tol );
            phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi; % + this.zero_tol);
            phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi; % + this.zero_tol);  

        end
        
        %% Returns the regressor vector 3rd time derivative ks*phi_3dot.
        %  @param[in] x: The phase variable.
        %  @param[in] x_dot: The phase variable 1st time derivative.
        %  @param[in] x_ddot: The phase variable 2nd time derivative.
        %  @param[in] x_3dot: The phase variable 3rd time derivative.
        %  @return regressor vector 3rd time derivative.
        function phi_3dot = regressVec3Dot(this, x, dx, ddx, d3x)
            
            if (x < this.x_min || x > this.x_max)
                phi_3dot = zeros(length(this.c),1);
                return;
            end
            
            psi = this.kernel_fun_ptr(x);
            psi_dot = this.kernelFunDot(x, dx);
            psi_ddot = this.kernelFunDDot(x, dx, ddx);
            psi_3dot = this.kernelFun3Dot(x, dx, ddx, d3x);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
            sum_psi_ddot = sum(psi_ddot);
            sum_psi_3dot = sum(psi_3dot);

            phi = psi / sum_psi;
            phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi;
            phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;
            phi_3dot = (psi_3dot - 3*phi_ddot*sum_psi_dot - 3*phi_dot*sum_psi_ddot - phi*sum_psi_3dot) / sum_psi;

        end

        %% Plots the kernels for the given input values.
        %  @param[in] x_data: vector of input values.
        function plotPsi(this, x_data)
            
            Psi = this.kernel_fun_ptr(x_data);
            figure;
            hold on;
            for i=1:length(this.c)
               plot(x_data, Psi(i,:), 'LineWidth',2); 
            end
            hold off

        end
        
        function plotRegressVec(this, x_data)
            
            n_ker = length(this.c);
            n_data = length(x_data);

            Phi = zeros(n_ker, n_data);
            Phi_dot = zeros(n_ker, n_data);
            Phi_ddot = zeros(n_ker, n_data);
            for j=1:n_data
                Phi(:,j) = this.regressVec(x_data(j));
                Phi_dot(:,j) = this.regressVecDot(x_data(j),1);
                Phi_ddot(:,j) = this.regressVecDDot(x_data(j),1,0);
            end
            
%             temp = sum(Phi,1);
%             figure;
%             plot(x_data, temp); % must be ones
            y_data = {Phi, Phi_dot, Phi_ddot};
            figure;
            for k=1:3
                subplot(3,1,k);
                hold on;
                for i=1:length(this.c)
                   plot(x_data, y_data{k}(i,:), 'LineWidth',2); 
                end
                axis tight;
                hold off;
            end

        end
        
    end
    
    %% =======================================================
    %% ===============  Protected functions  =================
    %% =======================================================
    methods (Access = protected)
        
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
        
        function psi = truncGaussKernel(this, x, zero_tol)

            psi = this.kernelFun(x);
            psi(psi<zero_tol) = 0;

        end
        
        function psi_dot = kernelFunDot(this, x, x_dot)

            n = length(x);
            psi = this.kernel_fun_ptr(x);
            psi_dot = zeros(length(this.c), n);
            
            for j=1:n
                a = (x(j)-this.c)*x_dot(j);
                psi_dot(:,j) = -2*this.h.*( psi(:,j).*a);
            end 

        end
        
        function psi_ddot = kernelFunDDot(this, x, x_dot, x_ddot)

            n = length(x);
            psi = this.kernel_fun_ptr(x);
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
            psi = this.kernel_fun_ptr(x);
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
    
    %% =====================================================
    %% ===============  Public properties  =================
    %% =====================================================
    
    properties (Access = public)

        c % N_kernels x 1 vector with the kernels' centers
        h % N_kernels x 1 vector with the kernels' inverse width
        
        kernel_std_scaling % the scaling of the kernel's std
        kernel_std % N_kernels x 1 vector with the kernels' std
        
    end
    
    %% ======================================================
    %% ===============  Private properties  =================
    %% ======================================================
    properties (Access = private)
        
        kernel_fun_ptr
        
        % range of input x outside which psi = 0 due to finite numerical precision
        x_min
        x_max 
        
    end
    
%     properties (Constant, Access = private)
%         
%         zero_tol = 1e-200 % small value used to avoid divisions with very small numbers
%         
%     end
    
end
