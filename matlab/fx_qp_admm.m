%% QP SOLVER - ADMM
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 16-02-2023
% Based on the work by Juan David Escárate
% ===============================================================================

function [tk, zk, uk] = fx_qp_admm(H, h, M, c, init_tk, init_zk, init_uk, opt_rho, opt_iters)          
%     global QP_vars counter
    % ------ QP PROBLEM ------ %
    % Min: 0.5*tk'*H*tk + h'*tk
    % ST:  M*tk <= c
    % 
    % ADMM Reformulation
    % Min: 0.5*tk'*H*tk + h'*tk + g(zk)
    % ST:  M*tk + zk = c
    % g(zk) is the indicator function of Z:
        % g(zk) = 0 if zk in Z,
        % g(zk) = ∞ if any component of zk not in Z
    %
    % Q -> R^(n x n)
    % A -> R^(n x n)
    % q -> R^(n)
    % c -> R^(n)
    % x -> R^(n)
    % z -> R^(n) (slack variable)
    
    % init_tk -> Initial tk point
    % init_zk -> Initial zk point
    % init_uk -> Initial uk point
    % rho -> Rho value for Augmented Lagrangian
    % iters -> Maximum iterations
    
    % ------ CHECK NUMBER OF INPUT ARGUMENTS ------ %
    if nargin < 4 || nargin > 9
        errorStruct.message = 'Incorrect number of inputs';
        errorStruct.identifier = 'ADMM_QP:incorrectInputs';
        error(errorStruct);
    end
    
    % Fill in unset optional values.
    %test_rho(A,Q);
    rho = dhang_rho(M,H);
%     QP_vars.rho(counter) = rho;
    iters = 50;
    switch nargin
        case 8
            rho = opt_rho;
        case 9
            rho = opt_rho;
            iters = opt_iters;
    end
    
    % ------ DATA INITIALIZATION ------ %
    % Length of variables
    n = size(h,1);
    % Quantity of restrictions
    n_c = size(c,1);
    % Data initialization
    if nargin > 4
        tk = init_tk;
        zk = init_zk;
        uk = init_uk;
    else
        tk = zeros(n,1);
        zk = zeros(n_c,1);
        uk = zeros(n_c,1);
    end
    % ------ ADMM FORMULATION ------ %

    % Augmented Lagrangian scaled form
    % L_rho(x, z, u) = f(x) + g(z) + (rho/2)*||Ax + z - c + u||_2^2
    % where u -> y/rho (scaled dual variable)

    % X minimization
        % To minimize we should calculate d(L_rho)/dx derivative
        % d(L_rho)/dx = f'(x) + rho*A^T(Ax + v) = 0
        % where v -> z - c + u
        % d(L_rho)/dx = Qx + q^T + rho*A^T(Ax + v) = 0
        % x = (Q + rho*A^T*A)^-1 * (-rho*A^T*v - q^T)

    % Z minimization
        %Usually evaluating Projection_Z(z) is inexpensive; 
        % for example, if Z = [α, β] is an interval,
        % Projection_Z(z) = min{max{z,α},β}
        % In this case Z = [0, inf) so
        % Projection_Z(z) = max{z,0}
        
    R = H + rho*(M'*M);
    R_inv = R \ eye(size(R,1));
%     if(counter == 1)
%         QP_vars.R_inv = R_inv;
%     end
    % Iterations
    for k = 1:iters
        v_x = zk - c + uk;
        tk = R_inv * (-rho*M'*v_x - h);    % update x
           
        zk = max(0, -M*tk - uk + c);             % update z

        % Then we update the scaled dual variable
        uk = uk + (M*tk + zk - c);                % update u
    end

end

function rho = dhang_rho(M,H)
    % ------ QP PROBLEM ------ %
    % Min: 1/2*x'*H*x + q'*x + g(z)
    % ST:  Mx + z = c
    %      z >= 0
    % g(z) is the indicator function of Z:
        % g(z) = 0 if z in Z,
        % g(z) = ∞ if any component of z not in Z

    % Define the singular decomposition of M as
    % M = U*S*V'
    % where U and V are orthonormal and S is diagonal with positive
    % real matrix entries
    [~,S,V] = svds(M);
    % Set Sd = (S'*S)^(-1/2)
    Sd = (S'*S)^(-1/2);
    % Calculate Pd as 
    % Pd = Sd'*V'*Q*V*Sd
    Pd = Sd'*V'*H*V*Sd;
    % Calculate non-zero eigenvalues
    lambda_Pd  = eig(Pd);
    lambda_max = max(lambda_Pd); 
    beta_1     = 1/lambda_Pd(1);
    beta_2     = 1/lambda_max;
    % Solve eq 36 from Dang paper (Embedded ADMM-based 
    % QP Solver for MPC with polytopic constraints)
    % a4*rho^4 + a3*rho^3 + a2*rho^2 + a1*rho + a0
    a4 =  beta_1*beta_2*(beta_1^2 + beta_2^2);
    a3 = -(beta_1^3 + beta_2^3 - 3*beta_1*beta_2*(beta_1 + beta_2));
    a2 = -(beta_1 - beta_2)^2;
    a1 = -2*(beta_1 + beta_2);
    a0 = -2;
    eq_roots = roots([a4, a3, a2, a1, a0]);
    rho = abs(eq_roots(1));
end