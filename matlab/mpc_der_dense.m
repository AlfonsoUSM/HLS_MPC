%% CONTROL MPC DE FILTRO LC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 06-10-2023
% Based on the work by Juan David Escárate
% ===============================================================================

clc; clear;

%% System

format('longE')

N_HOR = 2;      % tamaño del horizonte de predicción
x0   = [2; 0; 40; 0]; % Estado inicial

% parámetros del sistema
Rf  = 0.065;    % Resistencia Filtro 
Lf  = 3e-3;     % Inductancia Filtro 
Cf  = 15e-6;    % Capacitancia Filtro 
Vdc = 100;      % Voltaje DC de entrada, restricciones de voltaje en dq
Ts  = 200e-6;   % Periodo de muestreo
f   = 50;       % Frecuencia Hz
Wb  = 2*pi*f;   % Frecuencia rad/s   
RL  = 23.6;     % Resistencia de carga
Inom = 8;       % Corriente nominal, restricciones de corriente en dq
tsimu = 0.03;          % tsimu: Tiempo de simulación en segundos
k=0:Ts:tsimu-Ts;    % t: Arreglo de tiempo

% Modelo espacio estado filtro LC (Continuo)
%       x' = Ac*x + Bc*u + Bpc*d
%       y  = C*x
%
%       x = [id; iq; vd; vq]
%       d = [iod; ioq]
%       u = -k*x_nau + u_s (Actuación MPC)
Ac   = [-Rf/Lf      Wb  -1/Lf       0;
           -Wb  -Rf/Lf      0   -1/Lf;
          1/Cf       0      0      Wb;
             0     1/Cf   -Wb       0];

Bc   = [  1/Lf      0;      
             0   1/Lf;
             0      0;
             0      0];

Bpc  = [    0      0;      
            0      0;
        -1/Cf      0;
            0  -1/Cf];

C   = eye(4);

% Modelo espacio estado filtro LC (Discreto)
% x(k+1) = A*x(k) + B*u(k) + Bp*d(k) -> d = perturbacion (io, corriente de salida)
% y(k)   = Cx(k)
[~, Bp] = c2d(Ac,Bpc,Ts);
[A, B] = c2d(Ac,Bc,Ts);

%% Constraints

% Restricciones de corriente (politopo de 10 lados)
Hi   = [3.078 1; % Matriz de restrcciones que define el politopo
       -3.078 1; % Se puede restringir con 5 rectas
        0.726 1; 
       -0.726 1; 
            0 1];
Mi   = [1 0 0 0; % Matriz para selecionar id e iq
        0 1 0 0];
ang  = 36*pi/180; % pi/5
xmax = Inom*[               3.078;
                            3.078;
      (sin(ang) + 0.726*cos(ang)); 
      (sin(ang) + 0.726*cos(ang));
                       sin(2*ang)];
xmin = -xmax;

% Restricciones de voltaje (politopo de 10 lados)
Hv    = [3.078 1; % Matriz de restrcciones que define el politopo
        -3.078 1; % Se puede restringir con 5 rectas
         0.726 1; 
        -0.726 1; 
             0 1];

umax = Vdc/sqrt(3)*[            3.078;
                                3.078;
          (sin(ang) + 0.726*cos(ang)); 
          (sin(ang) + 0.726*cos(ang));
                          sin(2*ang)];
umin = -umax;

% Matrices de Peso 
Gamma      = 100*eye(2); %20*eye(2);    % 20 o 100 es por tuning
Omega      = eye(4);
Omega(1,1) = 1e2;       % Se le da mayor peso a los estados id iq
Omega(2,2) = 1e2;

% Solución a LQR discreto
[~,OmegaN,~] = dlqr(A,B,Omega,Gamma);

% r: salida deseada del sistema

ADMM_iters = 10;

%% Dense Formulation

N_SYS = size(A,1);      % numero de estados
M_SYS = size(B,2);      % numero de actuaciones
P_SYS = size(C,1);      % numero de salidas
C_SYS = 2*(size(xmin,1)+size(umin,1));  % numero de restricciones 

N_QP = N_HOR * M_SYS;
M_QP = N_HOR * C_SYS;
xk = zeros(N_SYS, length(k), 'single');
uk = zeros(M_SYS, length(k), 'single');
xk(:,1) = x0;

theta = zeros(N_QP, length(k), 'single');

[D,E] = fx_dense_matrices(A,B,N_HOR);       % constantes del sistema
K = blkdiag(kron(eye(N_HOR-1),Omega),OmegaN);
L = kron(eye(N_HOR),Gamma); 
Q = 2*(L+E'*K*E);                           % constante del sistema
Q = (Q+Q')/2;
Hi = kron(eye(N_HOR),Hi*Mi);
Hv = kron(eye(N_HOR),Hv);
F = [Hi*E;-Hi*E];
J = [Hv;-Hv];
G = single(2*D'*K*E);                       % constante del sistema
D = single(Hi*D);
a = single(kron(ones(N_HOR,1),umin));       % constante del sistema
b = single(kron(ones(N_HOR,1),umax));       % constante del sistema
d = single(kron(ones(N_HOR,1),xmin));       % constante del sistema
e = single(kron(ones(N_HOR,1),xmax));       % constante del sistema
H = [F;J];                                  % constante del sistema

rho = single(fx_dhang_rho(Q,H));
%rho = single(0.10070947);%single(62.963413);%

Q = single(Q);
H = single(H);

%% MPC Iteration

t_ADMM = zeros(N_QP, 1, 'single');
z_ADMM = zeros(M_QP, 1, 'single');
u_ADMM = zeros(M_QP, 1, 'single');

for i=1:length(k)
    q = (xk(:,i)'*G)';
    f = [e-D*xk(:,i); D*xk(:,i)-d];
    h = [f; b; -a];
    [t_ADMM, z_ADMM, u_ADMM] = fx_qp_admm(Q, q, H, h,t_ADMM, z_ADMM, u_ADMM, rho, ADMM_iters);
    uk(:,i) = t_ADMM(1:M_SYS);
    xk(:,i+1) = A*xk(:,i)+B*uk(:,i); % Cálculo del siguiente estado
    theta(:,i) = t_ADMM(:,1);
end

R = Q + rho*(H'*H);
R_inv = R \ eye(size(R,1));
W = -rho*H';            % RhoHt_neg

%% Plot

figure
plot(xk(1,:))
hold on
plot(xk(2,:))
plot(xk(3,:))
plot(xk(4,:))
plot(uk(1,:))
plot(uk(2,:))
grid on


%% Generate C++ file with Global Variables (constants)

txtfile = "samples/MPC_der_dense_N"+N_HOR+".cpp";
txtfileID = fopen(txtfile,'w');

fprintf(txtfileID, "\n#include "+char(34)+"system.hpp"+char(34)+"\n\n// HOR = "+N_HOR+"\n#if defined DENSE\n\n");

fx_cpp_print_matrix(txtfileID, H, "data_t H[M_QP][N_QP]", M_QP, N_QP)
fx_cpp_print_matrix(txtfileID, -a, "data_t a_neg[5*N_HOR]", 5*N_HOR)
fx_cpp_print_matrix(txtfileID, b, "data_t b[5*N_HOR]", 5*N_HOR)
fx_cpp_print_matrix(txtfileID, d, "data_t d[5*N_HOR]", (5*N_HOR))
fx_cpp_print_matrix(txtfileID, e, "data_t e[5*N_HOR]", (5*N_HOR))
fx_cpp_print_matrix(txtfileID, D, "data_t D[5*N_HOR][N_SYS]", (5*N_HOR), N_SYS)
fx_cpp_print_matrix(txtfileID, G, "data_t G[N_SYS][N_QP]", N_SYS, N_QP)

fx_cpp_print_matrix(txtfileID, R_inv, "data_t R_inv[N_QP][N_QP]", N_QP, N_QP)
fx_cpp_print_matrix(txtfileID, W, "data_t W[N_QP][M_QP]", N_QP, M_QP)

fprintf(txtfileID, "\n#else\n\n// SPARSE\n\n#endif");

fclose(txtfileID);


%% Generate .bin file with samples

binfile = "samples/MPC_der_dense_N"+N_HOR+".bin";
binfileID = fopen(binfile,'w');

nSamples = length(k);
data_t = 'single';
A = single(A);
B = single(B);

fwrite(binfileID, N_SYS,'uint8');
fwrite(binfileID, M_SYS,'uint8');
fwrite(binfileID, P_SYS,'uint8');
fwrite(binfileID, N_HOR,'uint8');
fwrite(binfileID, N_QP,'uint16');
fwrite(binfileID, M_QP,'uint16');
fwrite(binfileID, ADMM_iters,'uint16');
fwrite(binfileID, nSamples,'uint16'); 

% fwrite(binfileID,reshape(xmin',1,[]),data_t);
% fwrite(binfileID,reshape(xmax',1,[]),data_t);
% fwrite(binfileID,reshape(umin',1,[]),data_t);
% fwrite(binfileID,reshape(umax',1,[]),data_t);
% fwrite(binfileID,reshape(Q',1,[]),data_t);
% fwrite(binfileID,reshape(H',1,[]),data_t);
fwrite(binfileID,reshape(A',1,[]),data_t);
fwrite(binfileID,reshape(B',1,[]),data_t);
% fwrite(binfileID,rho,data_t);
% fwrite(binfileID,reshape(R_inv',1,[]),data_t);
% fwrite(binfileID,reshape(RhoHt_neg',1,[]),data_t);
 
for sample = 1:nSamples
%         fwrite(fileID,reshape(c_hat(:,sample)',1,[]),'double');
    fwrite(binfileID,reshape(xk(:,sample)',1,[]),data_t);
    fwrite(binfileID,reshape(uk(:,sample)',1,[]),data_t);
%         fwrite(fileID,reshape(theta(:,sample),1,[]),'double');
end

fclose(binfileID);


