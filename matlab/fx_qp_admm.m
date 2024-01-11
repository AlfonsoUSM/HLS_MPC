%% QP SOLVER - ADMM
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 16-02-2023
% Based on the work by Juan David Escárate
% ===============================================================================

function [tk, zk, uk] = fx_qp_admm(Q, q, G, g, init_tk, init_zk, init_uk, rho, iters)          
%     global QP_vars counter
% ------ QP PROBLEM ------ %
% Min: 0.5*tk'*Q*tk + q'*tk
% ST:  G*tk <= g
% 
% ADMM Reformulation
% Min: 0.5*tk'*Q*tk + q'*tk + h(zk)
% ST:  G*tk + zk = g
% h(zk) is the indicator function of Z:
%   h(zk) = 0 if zk in Z,
%   h(zk) = ∞ if any component of zk not in Z
%
% Q -> R^(n x n)
% G -> R^(n x n)
% q -> R^(n)
% g -> R^(n)
% x -> R^(n)
% z -> R^(n) (slack variable)

% init_tk -> Initial tk point
% init_zk -> Initial zk point
% init_uk -> Initial uk point
% rho -> Rho value for Augmented Lagrangian
% iters -> Maximum iterations
    
    % ------ CHECK NUMBER OF INPUT ARGUMENTS ------ %
    if nargin < 9 || nargin > 9
        errorStruct.message = 'Incorrect number of inputs';
        errorStruct.identifier = 'ADMM_QP:incorrectInputs';
        error(errorStruct);
    end

    % ------ DATA INITIALIZATION ------ %
    tk = init_tk;
    zk = init_zk;
    uk = init_uk;

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
        
    R = Q + rho*(G'*G);
    R_inv = R \ eye(size(R,1));
    % Iterations
    for k = 1:iters
        v_x = zk - g + uk;
        tk = R_inv * (-rho*G'*v_x - q);    % update x
        zk = max(0, -G*tk - uk + g);       % update z
        % Then we update the scaled dual variable
        uk = uk + (G*tk + zk - g);         % update u
    end

end