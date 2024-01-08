%% CONTROL MPC DE FILTRO LC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 07-01-2024
% Based on the work by Juan David Escárate
% ===============================================================================

clc; clear;

%% System

format('longE')

N_HOR = 2;      % tamaño del horizonte de predicción
x0   = [0; 0; 0; 0]; % Estado inicial
% rk: salida deseada del sistema
rk = repelem([0 50; 0 0], [1,1], [1,2500]);
% dk: perturbacion
load('ADMM_model_10_iters.mat', 'io')
dk = io';%*0.8;

% Parámetros del sistema
Rf  = 0.065;    % Resistencia Filtro 
Lf  = 3e-3;     % Inductancia Filtro 
Cf  = 15e-6;    % Capacitancia Filtro 
Vdc = 100;      % Voltaje DC de entrada, restricciones de voltaje en dq
Ts  = 200e-6;   % Periodo de muestreo
fb  = 50;       % Frecuencia Hz
Wb  = 2*pi*fb;   % Frecuencia rad/s   
RL  = 23.6;     % Resistencia de carga
Inom = 8;       % Corriente nominal, restricciones de corriente en dq
tsimu = 0.5;       % tsimu: Tiempo de simulación en segundos
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

Bdc  = [    0      0;      
            0      0;
        -1/Cf      0;
            0  -1/Cf];

C   = [0 0 1 0; 0 0 0 1];%eye(4);

% Modelo espacio estado filtro LC (Discreto)
% x(k+1) = A*x(k) + B*u(k) + Bp*d(k) -> d = perturbacion (io, corriente de salida)
% y(k)   = Cx(k)
[~, Bd] = c2d(Ac,Bdc,Ts);
[A, B] = c2d(Ac,Bc,Ts);

% Constraints
% Restricciones de las entradas, voltajes [Vmd; Vmq] (politopo de 10 lados)
H = [ 3.078 1; % Matriz de restrcciones que define el politopo
            -3.078 1; % Se puede restringir con 5 rectas
             0.726 1; 
            -0.726 1; 
             0 1];
ang  = 36*pi/180; % pi/5
umax = Vdc/sqrt(3)*[3.078;
                        3.078;
                        (sin(ang) + 0.726*cos(ang)); 
                        (sin(ang) + 0.726*cos(ang));
                        sin(2*ang)];
umin = -umax;

% Restricciones de estados, sólo corriente [Ifd; Ifq] (politopo de 10 lados)
I = [ 3.078 1 0 0; % Matriz de restrcciones que define el politopo
            -3.078 1 0 0; % Se puede restringir con 5 rectas
             0.726 1 0 0; % Matriz para selecionar id e iq
            -0.726 1 0 0; 
             0 1 0 0];
xmax = Inom*[   3.078;
                3.078;
                (sin(ang) + 0.726*cos(ang)); 
                (sin(ang) + 0.726*cos(ang));
                sin(2*ang)];
xmin = -xmax;

% Matrices de Peso 
Gamma      = 100*eye(2); %20*eye(2);    % 20 o 100 es por tuning
Omega      = eye(4);
Omega(1,1) = 1e2;       % Se le da mayor peso a los estados id iq
Omega(2,2) = 1e2;

% Solución a LQR discreto
[~,OmegaN,~] = dlqr(A,B,Omega,Gamma);

IT_ADMM = 10;

%% Dense Formulation

N_SYS = size(A,1);      % numero de estados
M_SYS = size(B,2);      % numero de actuaciones
P_SYS = size(C,1);      % numero de salidas
A_SYS = size(umax,1);   % numero de pares de restricciones de la entrada
B_SYS = size(xmax,1);   % numero de pares de restricciones del estado
C_SYS = 2*(A_SYS+B_SYS);% numero total de restricciones 

N_QP = N_HOR * M_SYS;
M_QP = N_HOR * C_SYS;
xk = zeros(N_SYS, length(k), 'single');
uk = zeros(M_SYS, length(k), 'single');
xk(:,1) = x0;

%theta = zeros(N_QP, length(k), 'single');

[D,E] = fx_dense_matrices(A,B,N_HOR);       % constantes del sistema
K = blkdiag(kron(eye(N_HOR-1),Omega),OmegaN);
L = kron(eye(N_HOR),Gamma); 
Q = 2*(L+E'*K*E);                           % constante del sistema
Q = (Q+Q')/2;
V = kron(eye(N_HOR),H);
W = kron(eye(N_HOR),I);

F = single(2*D'*K*E);                       % constante del sistema
D = single(D);
G = [V;-V;W*E;-W*E];                        % constante del sistema

rho = single(62.963413);%single(fx_dhang_rho(Q,G));%single(0.10070947);%

Q = single(Q);
G = single(G);


%% MPC Iteration

t_ADMM = zeros(N_QP, 1, 'single');
z_ADMM = zeros(M_QP, 1, 'single');
u_ADMM = zeros(M_QP, 1, 'single');

for i=1:length(k)
    [xinf, uinf] = fx_stationary(A, B, C, rk(:,i), Bd, dk(:,i));
    q = ((xk(:,i)-xinf)'*F)';
    Vu = single(V*kron(ones(N_HOR,1),uinf));
    WDx = single(W*(D*xk(:,i) + kron(ones(N_HOR,1),xinf)));
    c = kron(ones(N_HOR,1),single(umax));       % constante del sistema
    d = kron(ones(N_HOR,1),single(umin));       % constante del sistema
    e = kron(ones(N_HOR,1),single(xmax));       % constante del sistema
    f = kron(ones(N_HOR,1),single(xmin));       % constante del sistema
    g = [c-Vu; Vu-d; e-WDx; WDx-f];
    [t_ADMM, z_ADMM, u_ADMM] = fx_qp_admm(Q, q, G, g,t_ADMM, z_ADMM, u_ADMM, rho, IT_ADMM);
    uk(:,i) = t_ADMM(1:M_SYS)+uinf;
    xk(:,i+1) = A*xk(:,i)+B*uk(:,i)+Bd*dk(:,i); % Cálculo del siguiente estado
    %theta(:,i) = t_ADMM(:,1);
end

R = Q + rho*(G'*G);
R_inv = R \ eye(size(R,1));
P = -rho*G';            % RhoHt_neg
T = [A-eye(N_SYS),B;C,zeros(M_SYS,M_SYS)];
T_inv = T \ eye(N_SYS+M_SYS);

%% Plot

figure
plot(rk(1,:))
hold on
plot(rk(2,:))
plot(sqrt(xk(1,:).*xk(1,:)+xk(2,:).*xk(2,:)))
plot(xk(3,:))
plot(xk(4,:))
plot(uk(1,:))
plot(uk(2,:))
grid on
legend('Vrefd', 'Vrefq', '|Ifdq| (x0 and x1)', 'Vcd (x2)', 'Vcq (x3)','Input u0','Input u1')


%% Generate C++ file with Global Variables (constants)

txtfile = "samples2/MPC_der_dense_N"+N_HOR+".cpp";
txtfileID = fopen(txtfile,'w');

fprintf(txtfileID, "\n#include "+char(34)+"system.hpp"+char(34)+"\n\n// HOR = "+N_HOR+"\n#if defined DENSE\n\n");

fx_cpp_print_matrix(txtfileID, G, "data_t G[M_QP][N_QP]", M_QP, N_QP)
fx_cpp_print_matrix(txtfileID, c, "data_t c[A_SYS*N_HOR]", A_SYS*N_HOR)
fx_cpp_print_matrix(txtfileID, d, "data_t d[A_SYS*N_HOR]", A_SYS*N_HOR)
fx_cpp_print_matrix(txtfileID, e, "data_t f[B_SYS*N_HOR]", (B_SYS*N_HOR))
fx_cpp_print_matrix(txtfileID, f, "data_t e[B_SYS*N_HOR]", (B_SYS*N_HOR))
fx_cpp_print_matrix(txtfileID, D, "data_t D[N_SYS*N_HOR][N_SYS]", (N_SYS*N_HOR), N_SYS)
fx_cpp_print_matrix(txtfileID, F, "data_t F[N_SYS][N_QP]", N_SYS, N_QP)

fx_cpp_print_matrix(txtfileID, umin, "data_t umin[A_SYS]", A_SYS)
fx_cpp_print_matrix(txtfileID, umax, "data_t umax[A_SYS]", A_SYS)
fx_cpp_print_matrix(txtfileID, xmin, "data_t xmin[B_SYS]", B_SYS)
fx_cpp_print_matrix(txtfileID, xmax, "data_t xmax[B_SYS]", B_SYS)
fx_cpp_print_matrix(txtfileID, V, "data_t V[A_SYS*N_HOR][N_QP]", A_SYS*N_HOR, N_QP)
fx_cpp_print_matrix(txtfileID, W, "data_t W[B_SYS*N_HOR][N_QP]", B_SYS*N_HOR, N_QP)

fx_cpp_print_matrix(txtfileID, R_inv, "data_t R_inv[N_QP][N_QP]", N_QP, N_QP)
fx_cpp_print_matrix(txtfileID, T_inv, "data_t T_inv[N_SYS+M_SYS][N_SYS+M_SYS]", N_SYS+M_SYS, N_SYS+M_SYS)
fx_cpp_print_matrix(txtfileID, P, "data_t W[N_QP][M_QP]", N_QP, M_QP)


fprintf(txtfileID, "\n#else\n\n// SPARSE\n\n#endif");

fclose(txtfileID);


%% Generate .bin file with samples

binfile = "samples2/MPC_der_dense_N"+N_HOR+".bin";
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
fwrite(binfileID, IT_ADMM,'uint16');
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
    fwrite(binfileID,reshape(rk(:,sample)',1,[]),data_t);
    fwrite(binfileID,reshape(uk(:,sample)',1,[]),data_t);
%         fwrite(fileID,reshape(theta(:,sample),1,[]),'double');
end

fclose(binfileID);


