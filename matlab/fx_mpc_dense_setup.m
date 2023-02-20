%% MPC DENSE FORMULATION SETUP WITH INPUT AND STATE CONSTRAINTS 
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 09-02-2023
% Based on the work by Andrew Morrison https://github.com/morrisort/embeddedMPC/
% ===============================================================================
 
function [Acal,Ocal,Q,H,M_hat]=fx_mpc_dense_setup(A,B,OmegaN,Omega,Gamma,N_HOR)
    % ---Calculos que se tienen que hacer una sola vez --------------------
    [Acal,Ocal] = AOcal(A,B,N_HOR);
    R=kron(eye(N_HOR),Gamma);
    Q=blkdiag(kron(eye(N_HOR-1),Omega),OmegaN);
    %Gmm = mdiag(Gamma,N);                  
    %Omg = mdiag(Omega,OmegaN,N);
    % Matriz H del funcional de costo VN(x_0,vec{u})
    H=2*R+2*Ocal'*Q*Ocal;              
    H=(H+H')/2;
    % Matriz M de las restricciones lineales Mu<=c
    M_tilde=[Ocal;-Ocal];  
    % Convertir restricciones caja en restricciones de Desigualdad
    M_hat=[M_tilde;eye(N_HOR);-eye(N_HOR)];
end

% Función para generar matrices diagonales por bloques
function [M] = mdiag(varargin) 
    if nargin==2
        A=varargin{1}; N=varargin{2}; M=A;    
        for i=1:1:N-1; M=blkdiag(M,A); end
    end
    if nargin==3
        A=varargin{1}; B=varargin{2}; N=varargin{3};  M=A;
        for i=1:1:N-2; M=blkdiag(M,A); end
        M=blkdiag(M,B);
    end
end
