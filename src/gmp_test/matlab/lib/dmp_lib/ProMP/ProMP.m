%% ProMP class
%  Weighted sum of Guassians.
%

classdef ProMP < handle
    
    properties (Constant)
          
          %% enum TRAIN_METHOD
          LWR = 201;
          LS = 203;

    end
       
    properties
        N_kernels % number of kernels (basis functions)

        w % N_kernels x 1 vector with the kernels' weights
        c % N_kernels x 1 vector with the kernels' centers
        h % N_kernels x 1 vector with the kernels' inverse width

        zero_tol % small value used to avoid divisions with very small numbers
        
        can_clock_ptr % pointer to canonical clock

    end

    methods
        %% Weighted Sum of Gaussians constructor.
        %  @param[in] N_kernels: The number of kernels.
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std. (optional, default=1.0)
        %
        function this = ProMP(N_kernels, can_clock_ptr, kernel_std_scaling)

            if (nargin < 3), kernel_std_scaling = 1.0; end
            
            this.N_kernels = N_kernels;
            this.can_clock_ptr = can_clock_ptr;
            
            this.zero_tol = 1e-100; %realmin;

            this.w = zeros(this.N_kernels,2);
            this.c = ((1:this.N_kernels)-1)'/(this.N_kernels-1);
            this.h = 1./(kernel_std_scaling*(this.c(2:end)-this.c(1:end-1))).^2;
            this.h = [this.h; this.h(end)];
            
        end

        %% Returns the number of kernels.
        %  @return The number of kernels.
        %
        function n_ker = getNumOfKernels(this)
            
            n_ker = length(this.w);
            
        end
        
        
        %% Trains the ProMP.
        %  @param[in] train_method: The training method (see dmp_::TRAIN_METHOD enum).
        %  @param[in] Time: Row vector with the timestamps of the training data points.
        %  @param[in] Fd: Row vector with the desired values.
        %  @param[in] train_error: Optinal pointer to return the training error as norm(F-Fd)/n_data.
        %
        function [train_error, F, dF] = train(this, train_method, Time, Fd, dFd)
            
            n_data = length(Time);

            t_end = Time(end);
            this.setTau(t_end);
            
            s = ones(1, n_data);
            Phi = zeros(this.N_kernels, n_data);
            
            % Phi_dot = zeros(this.N_kernels, n_data); 
            
            for j=1:n_data
                x = this.phase(Time(j));
                Phi(:,j) = this.kernelFunction(x);
                % Phi_dot(:,j) = this.kernelFunctionDot(x);
            end
            
%             figure;
%             hold on
%             for i=1:size(Phi,1), plot(Time, Phi(i,:)); end
%             hold off
%             
%             figure;
%             hold on
%             for i=1:size(Phi_dot,1), plot(Time, Phi_dot(i,:)); end
%             hold off
            
      
            if (train_method == ProMP.LWR)
                this.w(:,1) = LWR(Phi, s, Fd, this.zero_tol);
                this.w(:,2) = LWR(Phi, s, dFd, this.zero_tol);
            elseif (train_method == ProMP.LS)
                this.w(:,1) = leastSquares(Phi, s, Fd, this.zero_tol);
                this.w(:,2) = leastSquares(Phi, s, dFd, this.zero_tol);
            else
                error('[ProMP::train]: Unsopported training method...');
            end

            if (nargout > 0)
                F = zeros(size(Fd));
                dF = zeros(size(dFd));
                for j=1:size(F,2)
                    x = this.phase(Time(j));
                    F(j) = this.output(x);
                    dF(j) = this.outputDot(x);
                end
                train_error = norm(F-Fd)/length(F) + norm(dF-dFd)/length(dF);
            end

        end

        
        %% Returns a column vector with the values of the kernel functions.
        %  @param[in] x: The phase variable.
        %  @return: Column vector with the values of the kernel functions.
        %
        function psi = kernelFunction(this, x)

            n = length(x);
            psi = zeros(this.N_kernels, n);

            for j=1:n
                psi(:,j) = exp(-this.h.*((x(j)-this.c).^2));
            end 

        end
        
        function psi_dot = kernelFunctionDot(this, x)

            dx = this.phaseDot(x);
            
            n = length(x);
            psi = this.kernelFunction(x);
            psi_dot = zeros(this.N_kernels, n);
            
            for j=1:n
                psi_dot(:,j) = -2*this.h.*((x(j)-this.c)).*psi(:,j) * dx(j);
            end 

        end
         
        
        %% Returns the normalized weighted sum of the Gaussians for the given phase variable (time instant).
        %  @param[in] x: The phase variable.
        %  @param[out] f: The normalized weighted sum of the Gaussians.
        %
        function f = output(this,x)
            
            Psi = this.kernelFunction(x);         
            f = dot(Psi,this.w(:,1)) / (sum(Psi) + this.zero_tol); % add 'zero_tol' to avoid numerical issues
            
        end
        
        function f_dot = outputDot(this, x)
            
            Psi = this.kernelFunction(x);         
            f_dot = dot(Psi,this.w(:,2)) / (sum(Psi) + this.zero_tol); % add 'zero_tol' to avoid numerical issues
            
        end
        
        function f_ddot = outputDDot(this, x)
            
            Psi = this.kernelFunction(x);
            Psi_dot = this.kernelFunctionDot(x);
            sum_Psi = sum(Psi);
            sum_Psi_dot = sum(Psi_dot);
            
            Phi = ( Psi_dot*sum_Psi - Psi*sum_Psi_dot ) / ( sum_Psi^2 + this.zero_tol);
            
            f_ddot = dot(Phi,this.w(:,2));
            
        end
        
        %% Returns the phase variable.
        %  @param[in] t: The time instant.
        %  @param[out] x: The phase variable for time 't'.
        function x = phase(this, t)
            
            x = this.can_clock_ptr.getPhase(t);

        end
        
        
        %% Returns the derivative of the phase variable.
        %  @param[in] x: The phase variable.
        %  @param[out] dx: The derivative of the phase variable.
        function dx = phaseDot(this, x)
            
            dx = this.can_clock_ptr.getPhaseDot(x);

        end
        
        %% Returns the time scale of the DMP.
        %  @param[out] tau: The time scale of the this.
        function tau = getTau(this)

            tau = this.can_clock_ptr.getTau();

        end
        
        
        %% Sets the time scale of the DMP.
        function setTau(this, tau)

            this.can_clock_ptr.setTau(tau);

        end
        
    end
end
