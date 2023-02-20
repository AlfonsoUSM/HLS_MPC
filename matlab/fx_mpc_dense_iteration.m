%% MPC DENSE ITERATION WITH INPUT AND STATE CONSTRAINTS 
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 09-02-2023
% Based on the work by Andrew Morrison https://github.com/morrisort/embeddedMPC/
% ===============================================================================

function [u_star,saveMat]=fx_mpc_dense_iteration(x,r,xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Q,H,M_hat,N,qpSol,linSol,iterPDIP,iterMINRES, tol, saveMat)
    % quadprog options
    options =  optimset('Display','off');

    % Change of variables
    [x_inf,u_inf] = stationaryStateValues(A,B,C,r);
    
    a_tilde=ones(N,1)*(umin-u_inf);                    
    b_tilde=ones(N,1)*(umax-u_inf);
    x_tilde=x-x_inf;
    h_tilde=(2*x_tilde'*Acal'*Q*Ocal)';
    Xmin=kron(ones(N,1),xmin-x_inf); 
    Xmax=kron(ones(N,1),xmax-x_inf); 
    c_tilde=[(Xmax-Acal*x_tilde);-(Xmin-Acal*x_tilde)];
    % Convertir restricciones caja en Desigualdad
    c_hat=[c_tilde;b_tilde;-a_tilde];
    % Solve QP problem
    if strcmp(qpSol,'pdip')
        [u_tilde_star,~,~,~,saveMat]=fx_qp_pdip(H,h_tilde,M_hat,c_hat,iterPDIP,iterMINRES,linSol,tol,saveMat);
    end
    if strcmp(qpSol,'quadprog')
        u_tilde_star=quadprog(H,h_tilde,M_hat,c_hat,[],[],[],[],[],options);
    end
    
    % Change of variables
    u_star=u_tilde_star(1)+u_inf;
    % ---------------------------------------------------------------------
    if saveMat.PDIP
        i = saveMat.i;
        saveMat.h_tilde(:,:,(i-1)+1) = [h_tilde];
        saveMat.c_hat(:,:,(i-1)+1) = [c_hat];
        saveMat.u_tilde_star(:,:,(i-1)+1) = [u_tilde_star];
    end
end



