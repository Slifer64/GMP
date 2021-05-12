%% math_ namespace
% 

classdef math_

    methods (Static, Access = public)
        
        
        %% ====================================
        %% =======    Quaternion    ===========
        %% ====================================
        
        %% Calculates the product of two quaternions.
        %  Each quaternion must be a 4x1 vector [w; v], where w is the 
        %  scalar part and v:3x1 the vector part.
        %  @param[in] Q1: First quaternion
        %  @param[in] Q2: Second quaternion
        %  @return Q12: The quaternion product Q1*Q2
        function Q12 = quatProd(Q1, Q2)

            Q1 = Q1(:);
            Q2 = Q2(:);

            w1 = Q1(1);
            v1 = Q1(2:4);

            w2 = Q2(1);
            v2 = Q2(2:4);

            Q12 = zeros(4,1);
            Q12(1) = w1*w2 - v1'*v2;
            Q12(2:4) = w1*v2 + w2*v1 + cross(v1,v2);
            % if Q12(1)<0, Q12=-Q12; end

        end
        
        
        %% Calculates the inverse of a quaternion.
        %  @param[in] Q: quaternion in the form [w; v], where w is the 
        %                the scalar part and v:3x1 the vector part.
        %  @return invQ: The quaternion inverse of Q.
        function invQ = quatInv(Q)
        
            invQ = zeros(size(Q));

            n = (Q'*Q);
            if (n < 1e-16)
                invQ = zeros(4,1);
                return;
            end

            invQ(1) = Q(1);
            invQ(2:4) = -Q(2:4);

            invQ = invQ / n;

        end
        
        
        %% Calculates the logarithmic map of a quaternion.
        %  @param[in] Q: quaternion in the form [w; v], where w is the 
        %                the scalar part and v:3x1 the vector part.
        %  @param[in] zero_tol: Zero tollerance to avoid divisions by 0 (optional default=1e-16).
        %  @return q_log: The quaternion logarithm of 'Q'.
        function q_log = quatLog(Q, zero_tol)

            if (nargin < 2), zero_tol = 1e-16; end
            
            w = Q(1);
            v = Q(2:4);
            norm_e = norm(v);

            if (norm_e > zero_tol)
                theta = 2*real(atan2(norm_e,w));
                q_log = theta*v/norm_e;
            else
                q_log = zeros(size(v));
            end

        end
        

        %% Calculates the quaternion exponential map of a 3D vector.
        %  @param[in] q_log: 3x1 vector.
        %  @param[in] zero_tol: Zero tollerance to avoid divisions by 0 (optional default=1e-16).
        %  @return Q: The quaternion exp of 'q_log'.
        function Q = quatExp(q_log, zero_tol)

            if (nargin < 2), zero_tol = 1e-16; end

            norm_q_log = norm(q_log);
            theta = norm_q_log;

            Q = zeros(4,1);

            if (norm_q_log > zero_tol)
            Q(1) = cos(theta/2);
            Q(2:4) = sin(theta/2)*q_log/norm_q_log;
            else
              Q = [1 0 0 0]';
            end

        end
        
        
        %% Calculates the quaternion difference between two quaternions.
        %  Each quaternion must be a 4x1 vector [w; v], where w is the 
        %  scalar part and v:3x1 the vector part.
        %  @param[in] Q1: First quaternion.
        %  @param[in] Q2: Second quaternion.
        %  @return Qdiff: The quaternion difference Q1*inv(Q2)
        function Qdiff = quatDiff(Q1, Q2)
        
            Qdiff = math_.quatProd(Q1, math_.quatInv(Q2));

        end


        %% Calculates the conjugate of a quaternion.
        %  @param[in] Q: quaternion in the form [w; v], where w is the 
        %                the scalar part and v:3x1 the vector part.
        %  @return Qconj: The conjugate quaternion of Q.
        function Qconj = quatConj(Q)

            Qconj = zeros(size(Q));

            Qconj(1) = Q(1);
            Qconj(2:4) = -Q(2:4);

        end
        
        
        %% Returns a quatenion as a 4x4 matrix.
        %  This matrix can be used to express the quaternion product as a
        %  matrix-vector product, i.e. quatProd(Q1,Q2) = quat2mat(Q1)*Q2
        %  @param[in] Q: quaternion in the form [w; v], where w is the 
        %                the scalar part and v:3x1 the vector part.
        %  @return qmat: The matrix expression of a quaternion.
        function qmat = quat2mat(Q)

            qmat = zeros(4,4);

            w = Q(1);
            v = Q(2:4);

            qmat(:,1) = Q;
            qmat(1,2:4) = -v;
            qmat(2:4,2:4) = w*eye(3,3) + math_.vector2ssMatrix(v);

        end

        
        %% Expresses a unit quaternions a pseudo-position 3D vector.
        %  Assuming 'Qo' as the origin, the pseudo-position is calculated 
        %  from the quaternion logarithm of the difference between Q and
        %  Qo. That is, qPos = log( Q * inv(Qo) )
        %  The original unit quaternion from the pseudo-position is given 
        %  by: Q = exp(qPos) * Qo
        %  @param[in] Q: 4 X N matrix where each column is a unit quaternio
        %  @param[in] Qo: 4 x 1 vector which is the unit quaternion origin (optional, default = [1 0 0 0]').
        %  @return qPos: 3 X N matrix the j-th column is the pseudo-position of the  
        %                j-th quaternion of Q, using as origin the j-th quaternion of Qo.
        function qPos = quat2qpos(Q, Qo)

            if (nargin < 2), Qo = [1 0 0 0]'; end
            
            n = size(Q,2);

            qPos = zeros(3, n);

            for i=1:size(qPos,2)
                qPos(:,i) = math_.quatLog(math_.quatProd(Q(:,i), math_.quatInv(Qo)));
            end

        end
        
        
        %% ==============================
        %% =======    Misc    ===========
        %% ==============================
        
        
        %% Returns the 4x4 skew-symmetric matrix of a 3D vector.
        %  @param[in] v: 3D vector.
        %  @return ssMat: 4x4 skew-symmetric matrix of 'v'.
        function ssMat = vector2ssMatrix(v)

            ssMat(1,1) = 0;
            ssMat(2,2) = 0;
            ssMat(3,3) = 0;
            ssMat(1,2) = -v(3);
            ssMat(2,1) = v(3);
            ssMat(1,3) = v(2);
            ssMat(3,1) = -v(2);
            ssMat(3,2) = v(1);
            ssMat(2,3) = -v(1);

        end
        
        
        %%  TODO
        function a = get5thOrderParams(x1, x2, y1, y2, y1_dot, y2_dot, y1_ddot, y2_ddot)

            a = zeros(6,1);

            b = [y1; y2; y1_dot; y2_dot;  y1_ddot;  y2_ddot];
            A = [ 1  x1    x1^2    x1^3     x1^4      x1^5
                  1  x2    x2^2    x2^3     x2^4      x2^5
                  0   1    2*x1   3*x1^2   4*x1^3    5*x1^4
                  0   1    2*x2   3*x2^2   4*x2^3    5*x2^4
                  0   0     2      6*x1    12*x1^2   20*x1^3
                  0   0     2      6*x2    12*x2^2   20*x2^3
                 ];

            a = A\b;

        end
        

    end

end
