%% CONTROL MPC DE FILTRO LC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 11-01-2024
% Based on the work by Juan David Escárate
% ===============================================================================

clc; clear;

%% System

format('longE')

N_HOR = 2;      % tamaño del horizonte de predicción
IT_ADMM = 10;
x0   = [0; 0; 0; 0]; % Estado inicial
% rk: salida deseada del sistema
rk = repelem([0 100; 0 0], [1,1], [1,2500]);
% RL Resistencia de carga
RL = repelem([47.2 10 30 20], [500 1000 500 501]);

% Parámetros del sistema
Rf  = 0.065;    % Resistencia Filtro 
Lf  = 3e-3;     % Inductancia Filtro 
Cf  = 15e-6;    % Capacitancia Filtro 
Vdc = 200;      % Voltaje DC de entrada, restricciones de voltaje en dq
Ts  = 200e-6;   % Periodo de muestreo
fb  = 50;       % Frecuencia Hz
Wb  = 2*pi*fb;  % Frecuencia rad/s   
Inom = 8;       % Corriente nominal, restricciones de corriente en dq
tsimu = 0.5;        % tsimu: Tiempo de simulación en segundos
k=0:Ts:tsimu-Ts;    % t: Arreglo de tiempo

% Modelo espacio estado filtro LC (Continuo)
%       x' = Ac*x + Bc*u + Bpc*d
%       y  = C*x
%
%       x = [id; iq; vd; vq]
%       d = [iod; ioq]
%       u = -k*x_nau + u_s (Actuación MPC)
A   = [-Rf/Lf      Wb  -1/Lf       0;
           -Wb  -Rf/Lf      0   -1/Lf;
          1/Cf       0      0      Wb;
             0     1/Cf   -Wb       0];

B   = [  1/Lf      0;      
             0   1/Lf;
             0      0;
             0      0];

C = eye(4); % matriz de salidas

Cu   = [0 0 1 0; 0 0 0 1]; % matriz de referencias

Bp  = [    0      0;      
            0      0;
        -1/Cf      0;
            0  -1/Cf];

% dk = [0 0 1/RL 0; 0 0 0 1/RL]*xk

% Modelo espacio estado filtro LC (Discreto)
% x(k+1) = A*x(k) + B*u(k) + Bp*d(k) -> d = perturbacion (io, corriente de salida)
% y(k)   = Cx(k)
[~, Bpd] = c2d(A,Bp,Ts);
[Ad, Bd] = c2d(A,B,Ts);

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
Omega(1,1) = 1e2;       % Se le da mayor peso a los estados id iq (corrientes)
Omega(2,2) = 1e2;

% Solución a LQR discreto
[~,OmegaN,~] = dlqr(Ad,Bd,Omega,Gamma);

%% Dense Formulation

N_SYS = size(Ad,1);         % numero de estados
M_SYS = size(Bd,2);         % numero de actuaciones
P_SYS = size(Cu,1);         % numero de referencias (salidas)
A_SYS = size(umax,1);       % numero de pares de restricciones de la entrada
B_SYS = size(xmax,1);       % numero de pares de restricciones del estado
C_SYS = 2*(A_SYS+B_SYS);    % numero total de restricciones 
D_SYS = size(Bpd,2);        % numero de perturbaciones

N_QP = N_HOR * M_SYS;
M_QP = N_HOR * C_SYS;
xk = zeros(N_SYS, length(k), 'single');
uk = zeros(M_SYS, length(k), 'single');
dk = zeros(M_SYS, length(k), 'single');
xk(:,1) = x0;
dk(:,1) = [0 0 1/RL(1) 0; 0 0 0 1/RL(1)]*x0;
J = zeros(length(k));
%theta = zeros(N_QP, length(k), 'single');

[D,E] = fx_dense_matrices(Ad,Bd,N_HOR);       % constantes del sistema
K = blkdiag(kron(eye(N_HOR-1),Omega),OmegaN);
L = kron(eye(N_HOR),Gamma); 
Q = 2*(L+E'*K*E);                           % constante del sistema
Q = (Q+Q')/2;
V = kron(eye(N_HOR),H);
W = kron(eye(N_HOR),I);

F = single(2*D'*K*E);                       % constante del sistema
D = single(D);
G = [V;-V;W*E;-W*E];                        % constante del sistema

rho = single(fx_dhang_rho(Q,G));%single(62.963413);%single(0.10070947);%

Q = single(Q);
G = single(G);


%% MPC Iteration

t_ADMM = zeros(N_QP, 1, 'single');
z_ADMM = zeros(M_QP, 1, 'single');
u_ADMM = zeros(M_QP, 1, 'single');

for i=1:length(k)
    dk(:,i) = [0 0 1/RL(i) 0; 0 0 0 1/RL(i)]*xk(:,i);
    [xinf, uinf] = fx_stationary(Ad, Bd, Cu, rk(:,i), Bpd, dk(:,i));
    q = ((xk(:,i)-xinf)'*F)';
    Huinf = H*uinf;
    Ixinf = I*xinf;
    WDx = single(W*D*(xk(:,i)-xinf));
    c = kron(ones(N_HOR,1),single(umax-Huinf));       % 
    d = kron(ones(N_HOR,1),single(Huinf-umin));       % 
    e = kron(ones(N_HOR,1),single(xmax-Ixinf));       % 
    f = kron(ones(N_HOR,1),single(Ixinf-xmin));       % 
    g = [c; d; e-WDx; f+WDx];
    [t_ADMM, z_ADMM, u_ADMM] = fx_qp_admm(Q, q, G, g,t_ADMM, z_ADMM, u_ADMM, rho, IT_ADMM);
    uk(:,i) = t_ADMM(1:M_SYS) + uinf;
    xk(:,i+1) = Ad*xk(:,i)+Bd*uk(:,i)+Bpd*dk(:,i); % Cálculo del siguiente estado
    %theta(:,i) = t_ADMM(:,1);
    J(i) = t_ADMM'*Q*t_ADMM/2 + q'*t_ADMM;
end

R = Q + rho*(G'*G);
R_inv = R \ eye(size(R,1));
P = -rho*G';            % RhoGt_neg
T = [eye(N_SYS)-Ad,-Bd;Cu,zeros(M_SYS,M_SYS)];
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

%% Plot 2

figure
plot(uk(2,:), uk(1,:))
hold on
plot(Vdc/(sqrt(3))*sin(0:pi/100:2*pi),Vdc/(sqrt(3))*cos(0:pi/100:2*pi))
grid on
temp = -150:50:150;
for i = 1:A_SYS
    y = (-H(i,1)*temp + umax(i))/H(i,2);
    plot(y,temp)
end
hold off
axis([-150 150 -150 150])
%legend('|Vdq|', 'Inom')

figure
plot(xk(2,:), xk(1,:))
hold on
plot(8*sin(0:pi/100:2*pi),8*cos(0:pi/100:2*pi))
grid on
hold off
legend('I', 'Inom')

%% Generate C++ file with Global Variables (constants)

txtfile = "samples2/MPC_der_dense_N"+N_HOR+".cpp";
txtfileID = fopen(txtfile,'w');

fprintf(txtfileID, "\n#include "+char(34)+"system.hpp"+char(34)+"\n\n// HOR = "+N_HOR+"\n#if defined DENSE\n\n");

fx_cpp_print_matrix(txtfileID, Bpd, "data_t Bpd[N_SYS][D_SYS]", N_SYS, D_SYS)
fx_cpp_print_matrix(txtfileID, G, "data_t G[M_QP][N_QP]", M_QP, N_QP)
%fx_cpp_print_matrix(txtfileID, c, "data_t c[A_SYS*N_HOR]", A_SYS*N_HOR)
%fx_cpp_print_matrix(txtfileID, d, "data_t d[A_SYS*N_HOR]", A_SYS*N_HOR)
%fx_cpp_print_matrix(txtfileID, e, "data_t f[B_SYS*N_HOR]", (B_SYS*N_HOR))
%fx_cpp_print_matrix(txtfileID, f, "data_t e[B_SYS*N_HOR]", (B_SYS*N_HOR))
%fx_cpp_print_matrix(txtfileID, D, "data_t D[N_SYS*N_HOR][N_SYS]", (N_SYS*N_HOR), N_SYS)
fx_cpp_print_matrix(txtfileID, F, "data_t F[N_SYS][N_QP]", N_SYS, N_QP)

fx_cpp_print_matrix(txtfileID, umin, "data_t umin[A_SYS]", A_SYS)
fx_cpp_print_matrix(txtfileID, umax, "data_t umax[A_SYS]", A_SYS)
fx_cpp_print_matrix(txtfileID, xmin, "data_t xmin[B_SYS]", B_SYS)
fx_cpp_print_matrix(txtfileID, xmax, "data_t xmax[B_SYS]", B_SYS)
fx_cpp_print_matrix(txtfileID, H, "data_t H[A_SYS][M_SYS]", A_SYS, M_SYS)
fx_cpp_print_matrix(txtfileID, I, "data_t I[B_SYS][N_SYS]", B_SYS, N_SYS)
%fx_cpp_print_matrix(txtfileID, V, "data_t V[A_SYS*N_HOR][N_QP]", A_SYS*N_HOR, N_QP)
%fx_cpp_print_matrix(txtfileID, W, "data_t W[B_SYS*N_HOR][N_SYS*N_HOR]", B_SYS*N_HOR, N_QP)
fx_cpp_print_matrix(txtfileID, W*D, "data_t WD[B_SYS*N_HOR][N_SYS]", B_SYS*N_HOR, N_QP)

fx_cpp_print_matrix(txtfileID, R_inv, "data_t R_inv[N_QP][N_QP]", N_QP, N_QP)
fx_cpp_print_matrix(txtfileID, T_inv, "data_t T_inv[N_SYS+M_SYS][N_SYS+M_SYS]", N_SYS+M_SYS, N_SYS+M_SYS)
fx_cpp_print_matrix(txtfileID, P, "data_t P[N_QP][M_QP]", N_QP, M_QP)

fprintf(txtfileID, "\n#else\n\n// SPARSE\n\n#endif");

fclose(txtfileID);


%% Generate .bin file with samples

binfile = "samples2/MPC_der_dense_N"+N_HOR+".bin";
binfileID = fopen(binfile,'w');

nSamples = length(k);
data_t = 'single';
Ad = single(Ad);
Bd = single(Bd);

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
fwrite(binfileID,reshape(Ad',1,[]),data_t);
fwrite(binfileID,reshape(Bd',1,[]),data_t);
% fwrite(binfileID,rho,data_t);
% fwrite(binfileID,reshape(R_inv',1,[]),data_t);
% fwrite(binfileID,reshape(RhoHt_neg',1,[]),data_t);
 
for sample = 1:nSamples
%         fwrite(fileID,reshape(c_hat(:,sample)',1,[]),'double');
    fwrite(binfileID,reshape(xk(:,sample)',1,[]),data_t);
    fwrite(binfileID,reshape(rk(:,sample)',1,[]),data_t);
    fwrite(binfileID,reshape(dk(:,sample)',1,[]),data_t);
    fwrite(binfileID,reshape(uk(:,sample)',1,[]),data_t);
%         fwrite(fileID,reshape(theta(:,sample),1,[]),'double');
end

fclose(binfileID);


