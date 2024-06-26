%% MPC for DC-DC motor, sparse formulation
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 06-09-2023
% Based on the work by Andrew Morrison and Juan Escárate
% https://github.com/morrisort/embeddedMPC/
% ===============================================================================
clear, clc

%% System

format('longE')

N_HOR = 4;      % Ns: tamaños del horizonte de predicción

% Arreglo de tiempo
Ts = 0.001;     % Ts: Periodo de muestreo en segundos
tsim = 3;       % tsim: Tiempo de simulación en segundos
k = 0:Ts:tsim-Ts;       % k: Arreglo de tiempo

% Datos del servomotor en tiempo discreto
kappa = 39.08/27.92;
tau = 1/27.92; 
A = [exp(-Ts/tau),0; tau*(1-exp(-Ts/tau)),1];
B = [kappa*(1-exp(-Ts/tau));Ts*kappa+tau*kappa*(exp(-Ts/tau)-1)];
C = [0,1];
x0 = single([3.0;-1.0]);    % x0: velocidad y posicion angular inicial


% Datos de las restricciones
% umin, umax en Volts:
umin = single(-3);
umax = single(3);
% xmin, xmax: rad/s, rad
xmin = single([-5;-2]);
xmax = single([5;2]);

Gamma = 0.1;
Omega = C'*C;
[Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);

% % r: salida deseada del sistema
% f = 0.2; %hz
% r = square(2*pi*f*t);

ADMM_iters = 10;

%% Spares Formulation

N_SYS = size(A,1);      % numero de estados
M_SYS = size(B,2);      % umero de actuaciones
P_SYS = size(C,1);      % numero de salidas

N_QP = (N_HOR * (N_SYS + M_SYS) + N_SYS);
M_QP = (4 * ( N_HOR + 1) * N_SYS + 2 * N_HOR * M_SYS);
xk = zeros(N_SYS, length(k), 'single');
uk = zeros(M_SYS, length(k), 'single');
xk(:,1) = x0;

theta = zeros(N_QP, length(k), 'single');
Q = single(blkdiag(kron(eye(N_HOR),Omega), OmegaN, kron(eye(N_HOR), Gamma)));   % constante del sistema
q = zeros(N_QP, 1, 'single');                                                   % constante del sistema

F = single([blkdiag(kron(eye(N_HOR),-eye(N_SYS)),eye(N_SYS))+kron([zeros(1,N_HOR+1); ...
       eye(N_HOR),zeros(N_HOR,1)],A),kron([zeros(1,N_HOR);eye(N_HOR)],B)]);
D = [eye(N_SYS);-eye(N_SYS)];
d = [xmax; -xmin];
DN = D;
dN = d;
E = [eye(M_SYS); -eye(M_SYS)];
e = [umax; -umin];
G = single(blkdiag(kron(eye(N_HOR),D), DN, kron(eye(N_HOR), E)));

f = zeros(N_SYS*(N_HOR+1), 1, 'single');
g = single([kron(ones(N_HOR,1),d);dN;kron(ones(N_HOR,1),e)]);                   % constante del sistema
H = [G; F; -F];                                                                 % constante del sistema

rho = single(fx_dhang_rho(Q,H));
%rho = single(0.10070947);%single(62.963413);%

%% MPC Iteration

t_ADMM = zeros(N_QP, 1, 'single');
z_ADMM = zeros(M_QP, 1, 'single');
u_ADMM = zeros(M_QP, 1, 'single');

for i=1:length(k)
    f(1:N_SYS,1) = -xk(:,i);
    h = [g; f; -f];
    [t_ADMM, z_ADMM, u_ADMM] = fx_qp_admm(Q, q, H, h,t_ADMM, z_ADMM, u_ADMM, rho, ADMM_iters);
    uk(:,i) = t_ADMM((N_HOR+1)*N_SYS+1:(N_HOR+1)*N_SYS+M_SYS,1);
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
plot(uk(1,:))
grid on

%% Generate C++ file with Global Variables (constants)

txtfile = "samples/MPC_motor_sparse_N"+N_HOR+".txt";
txtfileID = fopen(txtfile,'w');

fx_cpp_print_matrix(txtfileID, q, "data_t q[N_QP]", N_QP)
fx_cpp_print_matrix(txtfileID, H, "data_t H[M_QP][N_QP]", M_QP, N_QP)
fx_cpp_print_matrix(txtfileID, g, "data_t g[2*N_QP]", (2*N_QP))

fx_cpp_print_matrix(txtfileID, R_inv, "data_t R_inv[N_QP][N_QP]", N_QP, N_QP)
fx_cpp_print_matrix(txtfileID, W, "data_t W[N_QP][M_QP]", N_QP, M_QP)

fclose(txtfileID);

%% Generate .bin file with samples

binfile = "samples/MPC_motor_sparse_N"+N_HOR+".bin";
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
fwrite(binfileID, 0,'uint16');

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

%%
% 
% A = single(A);
% B = single(B);
% 
% outputMat = 'samples/MPC_motor_sparse_N'+ string(N_HOR) +'.mat';
% matObj = matfile(outputMat,'Writable',true);
% matObj.N_SYS = N_SYS;
% matObj.M_SYS = M_SYS;
% matObj.P_SYS = P_SYS;
% matObj.iter = ADMM_iters;
% matObj.N_QP = N_QP;
% matObj.M_QP = M_QP;
% 
% matObj.umin = umin;
% matObj.umax = umax;
% matObj.xmin = xmin;
% matObj.xmax = xmax;
% matObj.Q = Q;
% matObj.q = q;
% matObj.H = H;
% matObj.g = g;
% % matObj.theta = theta;
% matObj.uk = uk;
% matObj.xk = xk;
% matObj.A = A;
% matObj.B = B;
% matObj.rho = rho;
% matObj.R_inv = R_inv;
% matObj.RhoHt_neg = RhoHt_neg;
% writeMPCSamples(N_HOR, outputMat, 'motor', 'sparse', 'single');
