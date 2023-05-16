

Ns = 4;

Ts = 0.001;         % Ts: Periodo de muestreo en segundos
tsimu = 3;         % tsimu: Tiempo de simulación en segundos
t=0:Ts:tsimu-Ts;    % t: Arreglo de tiempo

% Datos del servomotor en tiempo discreto
K=39.08/27.92;
tau=1/27.92; 
A = [exp(-Ts/tau),0; tau*(1-exp(-Ts/tau)),1];
B = [K*(1-exp(-Ts/tau));Ts*K+tau*K*(exp(-Ts/tau)-1)];
C = [0,1];
x0= [3;-1];          % x0: velocidad y posicion angular inicial

N_SYS = size(A,1);      % numero de estados
M_SYS = size(B,2);      % umero de actuaciones
P_SYS = size(C,1);      % numero de salidas

% Datos de las restricciones
% umin, umax en Volts:
umin=-10;
umax=10;
% xmin, xmax: rad/s, rad
xmin=[-5;-2];
xmax=[5;2];

Gamma=0.1;
Omega=C'*C;
[Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);

% % r: salida deseada del sistema
% f = 0.2; %hz
% r = square(2*pi*f*t);

ADMM_iters = 10;
rho = 62.963413;

%%
for N_HOR = Ns
    disp(['Procesando horizonte de tamaño: ', num2str(N_HOR)])
    N_QP = (N_HOR * (N_SYS + M_SYS) + N_SYS);
    M_QP = (4 * ( N_HOR + 1) * N_SYS + 2 * N_HOR * M_SYS);
    xk = zeros(N_SYS,length(t));
    uk = zeros(M_SYS,length(t));
    xk(:,1) = x0;

    % formulación sparse
    theta = zeros(N_QP, length(t));
%     c_hat = zeros(M_QP, length(t));
    Hcal = blkdiag(kron(speye(N_HOR),Omega), OmegaN, kron(speye(N_HOR), Gamma));
    h = zeros(N_QP,1);

    Fcal = [blkdiag(kron(speye(N_HOR),-eye(N_SYS)),eye(N_SYS))+kron([zeros(1,N_HOR+1);speye(N_HOR),zeros(N_HOR,1)],A),kron([zeros(1,N_HOR);speye(N_HOR)],B)];
    f = zeros(N_SYS*(N_HOR+1),1);
    J = [eye(N_SYS);-eye(N_SYS)];
    d = [xmax; -xmin];
    JN = J;
    dN = d;
    E = [eye(M_SYS);-eye(M_SYS)];
    e = [umax; -umin];
    Gcal = blkdiag(kron(speye(N_HOR),J), JN, kron(speye(N_HOR), E));
    g = [kron(ones(N_HOR,1),d);dN;kron(ones(N_HOR,1),e)];
    
    M_hat = [Gcal; Fcal; -Fcal];

    t_AMMD = zeros(N_QP,1);
    z_AMMD = zeros(M_QP,1);
    u_AMMD = zeros(M_QP,1);

    % setup

    for i=1:length(t)
        f(1:N_SYS,1) = -xk(:,i);
        c_hat = [g; f; -f];
        [t_AMMD, z_AMMD, u_AMMD] = fx_qp_admm(Hcal, h, M_hat, c_hat,t_AMMD, z_AMMD, u_AMMD, rho, ADMM_iters);
        uk(:,i) = t_AMMD((N_HOR+1)*N_SYS+1:(N_HOR+1)*N_SYS+M_SYS,1);
        xk(:,i+1) = A*xk(:,i)+B*uk(:,i); % Cálculo del siguiente estado
        theta(:,i) = t_AMMD(:,1);
    end

end

%%
figure
plot(xk(1,:))
hold on
plot(xk(2,:))
plot(uk(1,:))
grid on


%%



outputMat = 'samples/motor_MPC_N'+ string(N_HOR) +'.mat';
matObj = matfile(outputMat,'Writable',true);
matObj.N_SYS = N_SYS;
matObj.M_SYS = M_SYS;
matObj.P_SYS = P_SYS;
matObj.iter = ADMM_iters;
matObj.N_QP = N_QP;
matObj.M_QP = M_QP;

matObj.umin = umin;
matObj.umax = umax;
matObj.xmin = xmin;
matObj.xmax = xmax;
matObj.H = Hcal;
matObj.h = h;
matObj.M_hat = M_hat;
matObj.g = g;
% matObj.theta = theta;
matObj.uk = uk;
matObj.xk = xk;
matObj.A = A;
matObj.B = B;
writeMPCSamples(N_HOR, outputMat, 'motor', 'double');
