

Ns = 4;

Ts = 0.001;         % Ts: Periodo de muestreo en segundos
tsimu = 10;         % tsimu: Tiempo de simulación en segundos
t=0:Ts:tsimu-Ts;    % t: Arreglo de tiempo

% Datos del servomotor en tiempo discreto
K=39.08/27.92;
tau=1/27.92; 
A = [exp(-Ts/tau),0; tau*(1-exp(-Ts/tau)),1];
B = [K*(1-exp(-Ts/tau));Ts*K+tau*K*(exp(-Ts/tau)-1)];
C = [0,1];
x0= [0.3;0.1];          % x0: velocidad y posicion angular inicial

n = size(A,1);      % numero de estados
m = size(B,2);      % umero de actuaciones
p = size(C,1);      % numero de salidas

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

ADMM_iters = 25;

%%
for N_HOR = Ns
    disp(['Procesando horizonte de tamaño: ', num2str(N_HOR)])
    xk = zeros(n,length(t));
    uk = zeros(m,length(t));
    xk(:,1) = x0;

    % formulación sparse
    theta = zeros(N_HOR*(n+m)+n,length(t));
    Hcal = blkdiag(kron(speye(N_HOR),Omega), OmegaN, kron(speye(N_HOR), Gamma));
    h = zeros(N_HOR*(n+m)+n,1);

    Fcal = [blkdiag(kron(speye(N_HOR),-eye(n)),eye(n))+kron([zeros(1,N_HOR+1);speye(N_HOR),zeros(N_HOR,1)],A),kron([zeros(1,N_HOR);speye(N_HOR)],B)];
    f = zeros(n*(N_HOR+1),1);
    J = [eye(n);-eye(n)];
    d = [xmax; -xmin];
    JN = J;
    dN = d;
    E = [eye(m);-eye(m)];
    e = [umax; -umin];
    Gcal = blkdiag(kron(speye(N_HOR),J), JN, kron(speye(N_HOR), E));
    g = [kron(ones(N_HOR,1),d);dN;kron(ones(N_HOR,1),e)];
    
    M_hat = [Gcal; Fcal; -Fcal];

    theta_AMMD = zeros(N_HOR*(n+m)+n,1);
    z_AMMD = zeros(size(M_hat,1),1);
    u_AMMD = zeros(size(M_hat,1),1);

    % setup

    for i=1:length(t)
        f(1:n,1) = -xk(:,i);
        c_hat = [g; f; -f];
        [theta_AMMD, z_AMMD, u_AMMD] = fx_qp_admm(Hcal, h, M_hat, c_hat,theta_AMMD, z_AMMD, u_AMMD, 62.963413, ADMM_iters);
        uk(:,i) = theta_AMMD((N_HOR+1)*n+1:(N_HOR+1)*n+m,1);
        xk(:,i+1) = A*xk(:,i)+B*uk(:,i); % Cálculo del siguiente estado
        theta(:,i) = theta_AMMD(:,1);
    end

end

%%
figure
plot(xk(1,:))
hold on
plot(xk(2,:))
plot(uk(1,:))
grid on
