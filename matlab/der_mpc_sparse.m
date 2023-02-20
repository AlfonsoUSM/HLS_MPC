%% CONTROL MPC DE FILTRO LC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 10-02-2023
% Based on the work by Juan David Escárate
% ===============================================================================

%function ctrl = Sparse_Oficial_Lab(Input_funtion)
global x_ADMM z_ADMM u_ADMM QP_vars counter solver ADMM_iters
%% VALORES NOMINALES del Filtro:
% Valores de Planta laboratorio filtro LC
Rf  = 0.065;    % Resistencia Filtro 
Lf  = 3e-3;     % Inductancia Filtro 
Cf  = 15e-6;    % Capacitancia Filtro 
Vdc = 100;      % Voltaje DC de entrada
T   = 200e-6;   % Periodo de muestreo
f   = 50;       % Frecuencia Hz
Wb  = 2*pi*f;   % Frecuencia rad/s   
RL  = 23.6;     % Resistencia de carga

solver = "ADMM";
ADMM_iters = 10;
x_ADMM = zeros(16,1);
z_ADMM = zeros(74,1);
u_ADMM = zeros(74,1);
counter = 1;
max_steps = 2501;
QP_vars = {};
QP_vars.io = zeros(max_steps,2,1);
QP_vars.x0 = zeros(max_steps,4,1);
QP_vars.Q = zeros(16,16);
QP_vars.q = zeros(16,1);
QP_vars.A = zeros(74,16);
QP_vars.c = zeros(max_steps,74,1);
QP_vars.x_QP  = zeros(max_steps,16,1);
QP_vars.rho   = zeros(max_steps,1);

% Horizonte de PREDICCION
N = 2;

% Restricciones de voltaje y corriente en dq
Inom = 8; % Corriente nominal

% Modelo espacio estado filtro LC (Continuo)
% x' = A*x + B*u + Bp*d
% y  = C*x

% x = [id; iq; vd; vq]
% d = [iod; ioq]
% u = -k*x_nau + u_s (Actuación MPC)
A   = [ -Rf/Lf      Wb  -1/Lf       0;
           -Wb  -Rf/Lf      0   -1/Lf;
          1/Cf       0      0      Wb;
             0     1/Cf   -Wb       0];

B   = [   1/Lf      0;      
             0   1/Lf;
             0      0;
             0      0];

Bp  = [      0      0;      
             0      0;
         -1/Cf      0;
             0  -1/Cf];

C   = eye(4);

% Modelo espacio estado filtro LC (Discreto)
% x(k+1) = Ad*x(k) + Bd*u(k) + Bpd*d(k) -> d = perturbacion (io, corriente de salida)
% y(k)   = Cx(k)
[~, Bpd] = c2d(A,Bp,T);
[Ad, Bd] = c2d(A,B,T);
[nx, nu] = size(Bd);
n = length(Ad);

% Matrices de Peso 
%Gamma      = 20*eye(2); % 20 es por tuning
Gamma      = 100*eye(2); % 100 es por tuning
Omega      = eye(4);
% Se le da mayor peso a los estados id iq
Omega(1,1) = 1e2;
Omega(2,2) = 1e2;

% Solución a LQR discreto
[~,Omega_N,~] = dlqr(Ad,Bd,Omega,Gamma);

% Restricciones de corriente (politopo de 10 lados)
Hi   = [3.078 1; % Matriz de restrcciones que define el politopo
       -3.078 1; % Se puede restringir con 5 rectas
        0.726 1; 
       -0.726 1; 
            0 1];
Mi   = [1 0 0 0; % Matriz para selecionar id e iq
        0 1 0 0];
ang  = 36*pi/180; %
xmax = Inom*[               3.078;
                            3.078;
      (sin(ang) + 0.726*cos(ang)); 
      (sin(ang) + 0.726*cos(ang));
                       sin(2*ang)];
xmin = -xmax;

% Restricciones de voltaje (hexágono)
% Hv   = [       0  1; 
%     -3/sqrt(3) 1; 
%      3/sqrt(3) 1];
% 
% umax = Vdc*[ 1/sqrt(3);
%           2/sqrt(3); 
%           2/sqrt(3)];
% umin = -umax;
% Nueva restriccion de voltaje
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

% Cambio de nombre para matrices de peso
Q_N = Omega_N; % Region final
Q   = Omega;   % Estado 
R   = Gamma;   % Control

% Demultiplexación de entradas a la función
Vsd   = 50; %Input_funtion(1);
Vsq   = 0; %Input_funtion(2);
x_0   = [2; 0; 40; 0]; %Input_funtion(3:6);
i_od  = 0; %Input_funtion(7);
i_oq  = 0; %Input_funtion(8);
QP_vars.io(counter,:,:)    = [i_od;i_oq];
QP_vars.x0(counter,:,:)    = x_0;
%theta = Input_funtion(9);

% Referencia de salida
yref = [Vsd;  Vsq];
% Vd -> step; Vq -> 0

% Perturbaciones
d3   = 0;
d4   = 0;
per  = [i_od;
        i_oq;
          d3;
          d4];

% Obtener "x" y "u" en estado estacionario
[xinf, uinf] = Infinity(Ad, Bd, Bpd, C, per, yref);

%% Solver QP

% x = 0 -> N
% u = 0 -> N-1
% Con Horizonte = 2
% x_QP' -> [ x_0' x_1' x_2' u_0' u_1']
%       -> [ id_0 iq_0 vd_0 vq_0  id_1 iq_1 vd_1 vq_1  id_2 iq_2 vd_2 vq_2  Vsd_0 Vsq_0  Vsd_1 Vsq_1]

% Sólo existen restricciones en amplitud de Vsd, Vsq (Voltaje de entrada al 
% filtro) y id, iq (corrientes del filtro)

% Matrices Q y A de modelo QP
Q = blkdiag(kron(speye(N), Q), Q_N, kron(speye(N), R));
q = zeros(length(Q),1);

% Matriz de rotación
% T = [cos(theta) -sin(theta); 
%      sin(theta)  cos(theta)];

% Restricciones de igualdadd asociadas al modelo dinámico
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
Bu = kron([sparse(1, N); speye(N)], Bd);
% Se reescriben las restricciones de igualdad como desigualdad
% leq <= Aeq * x_QP <= ueq
Aeq = [Ax, Bu];
leq = [-(x_0 - xinf); zeros(N*nx, 1)];
ueq = leq;

% Restricciones de desigualdad asociada a las entradas y estados del
% sistema
% Aineq = blkdiag(kron(eye(N+1),Hi*Mi),kron(eye(N),Hv*T));
% uineq = [kron(ones(N+1,1),xmax-Hi*Mi*xinf);kron(ones(N,1),umax-Hv*T*uinf)];
% lineq = [kron(ones(N+1,1),xmin-Hi*Mi*xinf);kron(ones(N,1),umin-Hv*T*uinf)];

Aineq = blkdiag(kron(eye(N+1),Hi*Mi),kron(eye(N),Hv));
uineq = [kron(ones(N+1,1),xmax-Hi*Mi*xinf);kron(ones(N,1),umax-Hv*uinf)];
lineq = [kron(ones(N+1,1),xmin-Hi*Mi*xinf);kron(ones(N,1),umin-Hv*uinf)];

% Se eliminan las restricciones para los estados iniciales
Aineq(1:5,1:4) = zeros(5,4);
uineq(1:5)     = zeros(5,1);
lineq(1:5)     = zeros(5,1);

A  = [Aeq; Aineq];
lb = [leq; lineq];
ub = [ueq; uineq];

A  = [A;-A;];
c  = [ub;-lb];

if strcmp(solver, "quadprog")
    %%%%%%%%%%%%%%%%%%%% Quadprog:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = optimoptions('quadprog','Algorithm','interior-point-convex','LinearSolver','auto','Display','off');
    [x_QP]  = quadprog(Q,q,A,c,[],[],[],[],[],options);  
    ctrl    = x_QP((N+1)*nx+1:(N+1)*nx+nu)+uinf;

elseif strcmp(solver,"ADMM")
    %%%%%%%%%%%%%%%%%%%% ADMM:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x_ADMM,z_ADMM,u_ADMM]=fx_qp_admm(Q,q,A,c,x_ADMM,z_ADMM,u_ADMM,62.963413, ADMM_iters);
    ctrl = x_ADMM((N+1)*nx+1:(N+1)*nx+nu)+uinf;
    x_QP = x_ADMM;
else
    errorStruct.message = 'No solver defined in variable "solver"';
    errorStruct.identifier = 'Sparse_Oficial_Lab:nosolver';
    error(errorStruct);
end

%% %%%%%%%%%%%%%%%%%% Guardar Matrices y resultado:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if counter == 1
    QP_vars.Q = Q;
    QP_vars.q = q;
    QP_vars.A = A;
end
QP_vars.c(counter,:,:)    =    c;
QP_vars.x_QP(counter,:,:) = x_QP;
counter = counter + 1;

%end