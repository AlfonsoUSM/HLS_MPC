%% CONTROL MPC DE SERVOMOTOR CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 09-02-2023
% Based on the work by Andrew Morrison https://github.com/morrisort/embeddedMPC/
% ===============================================================================

% CONTROL MPC CON RESTRICCIONES EN LA ENTRADA Y EN LOS ESTADOS 
%   + SEGUIMIENTO SE SETPOINT
%clc; clear all; close all;

% Ns: tamaños del horizonte de predicción
% ejemplo: Ns = [2, 3, 4, 6, 8, 16, 32, 64];
Ns = 4;%[2, 3, 4, 6, 8, 16, 32, 64];

% Numero de iteraciones para PDIP
iterPDIP = 20;
% Linear solver: 'minres', 'cgrad', 'chol' o 'matlab'
linSol = "chol";
% Número de iteraciones para minres
iterMINRES = 30;
% tolerancia para minres y cgrad
tol = 1e-9;

% Comparar resultados con quadprog y LQR
compare = false;

% Guardar matrices de interés 
saveMat.MPC = false;% true;
saveMat.PDIP = false;% true;
saveMat.LS = false;% true; % esta opción hace que sea mucho mas lento el código

% la función controlMpc necesita saber en que iteración va para
% guardar las matrices en el lugar que corresponde
saveMat.i = 0;

% Donde guardar las matrices que van a ser utilizadas como referencia
saveMat.outputMat = 'data/servoMats_N';

% Arreglo de tiempo
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

% Datos de las restricciones
% umin, umax en Volts:
umin=-10;
umax=10;
% xmin, xmax: rad/s, rad
xmin=[-5;-2];
xmax=[5;2];

% Gamma y Omega se utilizan para MPC y LQR
Gamma=0.1;
Omega=C'*C;
[Linf,OmegaN,~]=dlqr(A,B,Omega,Gamma);

% r: salida deseada del sistema
f = 0.2; %hz
%r = sin(2*pi*f*t);
%r = -1 + (2)*rand(1,length(t));
r = 0*t;%square(2*pi*f*t);
%r=[ones(1,length(t)/4),-ones(1,length(t)/4), ...
%      ones(1,length(t)/4),-ones(1,length(t)/4)];


% Cálculo de la señal de control óptima y estado - control MPC
for N_HOR=Ns
    disp(['Procesando horizonte de tamaño: ', num2str(N_HOR)])
    
    % Reservar espacio memoria para guardar matrices ----------------------
    if saveMat.PDIP
        saveMat.h_tilde = zeros(N_HOR,1,length(t));
        saveMat.c_hat = zeros(6*N_HOR,1,length(t));
        saveMat.u_tilde_star = zeros(N_HOR,1,length(t));
    end
    if saveMat.LS
        %saveMat.Ak = zeros(n,n,length(t)*iterPDIP);
        %saveMat.bk = zeros(n,1,length(t)*iterPDIP);
        %saveMat.zk = zeros(n,1,length(t)*iterPDIP);
        
        saveMat.Ak = [];
        saveMat.bk = [];
        saveMat.zk = [];
    end
    
    
    % Control MPC ---------------------------------------------------------
    % Reservar espacio de memoria
    xPDIP = zeros(2,length(t));
    uPDIP = zeros(1,length(t));
    % Estado inicial
    xPDIP(:,1) = x0;

    % Preparar MPC
    [Acal,Ocal,Q,H,M_hat] = fx_mpc_dense_setup(A,B,OmegaN,Omega,Gamma,N_HOR);

%     profile on
    for i=1:length(t)
        saveMat.i = i;
        [uPDIP(i),saveMat]=fx_mpc_dense_iteration(xPDIP(:,i),r(i),xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Q,H,M_hat,N_HOR,'pdip',linSol,iterPDIP,iterMINRES, tol, saveMat);
        % Cálculo del siguiente estado
        xPDIP(:,i+1)=A*xPDIP(:,i)+B*uPDIP(1,i);
    end
    
%     p = profile("info");
%     profsave(p,"profiling/N_"+ string(N_HOR))
    
    xPDIP = xPDIP(:,1:end-1);

    
    % Guardar matrices de interés -----------------------------------------
    if saveMat.MPC || saveMat.PDIP || saveMat.LS
        % Si es que no existe el .mat, se crea
        outputMat = saveMat.outputMat + string(N_HOR) + '.mat';
        if exist(outputMat, 'file')
            save(outputMat,'N_HOR','-append');
        else
            save(outputMat,'N_HOR', '-v7.3');
        end
        matObj = matfile(outputMat,'Writable',true);
    end
    if saveMat.LS
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES;
        matObj.tol = tol;
        matObj.Ak = saveMat.Ak;
        matObj.bk = saveMat.bk;
        matObj.zk = saveMat.zk;
        writeLSSamples(N_HOR, outputMat);
    end
    if saveMat.PDIP
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES;
        matObj.iterPDIP = iterPDIP;
        matObj.tol = tol;
        matObj.H = H;
        matObj.M_hat = M_hat;
        matObj.x = xPDIP;
        matObj.h_tilde = saveMat.h_tilde;
        matObj.c_hat = saveMat.c_hat;
        matObj.u_tilde_star = saveMat.u_tilde_star;
        writePDIPSamples(N_HOR, outputMat);
    end
    if saveMat.MPC
        matObj.linSol = linSol;
        matObj.iterMINRES = iterMINRES;
        matObj.iterPDIP = iterPDIP;
        matObj.tol = tol;
        matObj.umin = umin;
        matObj.umax = umax;
        matObj.xmin = xmin;
        matObj.xmax = xmax;
        matObj.Acal = Acal;
        L=[eye(size(A,1))-A,-B;C,zeros(1,1)]; 
        L_inv = inv(L);
        matObj.L_invLast = L_inv(:,3);
        matObj.AcalQOcal = Acal'*Q*Ocal;
        matObj.Acal = Acal;
        matObj.r = r;
        matObj.H = H;
        matObj.M_hat = M_hat;
        matObj.u_star = uPDIP;
        matObj.x = xPDIP;
        matObj.A = A;
        matObj.B = B;
        writeMPCSamples(N_HOR, outputMat);
    end


    
    % Comparar resultados con Quadprog y LQR ------------------------------
    if compare
        % Reservar espacio de memoria
        xQuadProg = zeros(2,length(t));
        uQuadProg = zeros(N_HOR,length(t));

        % Estado inicial
        xQuadProg(:,1) = x0;

        saveMat.MPC = false;
        saveMat.PDIP = false;
        saveMat.LS = false;
        
        for i=1:length(t)
            [uQuadProg(:,i),saveMat]=fx_mpc_dense_iteration(xQuadProg(:,i),r(i),xmin,xmax,umin,umax,A,B,C,Acal,Ocal,Q,H,M_hat,N_HOR,'quadprog',linSol,iterPDIP,iterMINRES,tol,saveMat);
            % Cálculo del siguiente estado
            xQuadProg(:,i+1)=A*xPDIP(:,i)+B*uQuadProg(1,i);
        end
        xQuadProg = xQuadProg(:,1:end-1);
        % Evolución del sistema ante la entrada óptima obtenida con LQR:
        xLQR(:,1)=x0;
        L=[eye(size(A,1))-A,-B;C,zeros(1,1)]; 
        for i=1:length(t)
            bl=[zeros(size(A,1),1);r(1,i)];
            infy=L\bl; 
            xinfy=infy(1:size(A,1));         
            uinfy=infy(end);
            uLQR(1,i)=-Linf*(xLQR(:,i)-xinfy)+uinfy;
            xLQR(:,i+1) = A*xLQR(:,i) + B*uLQR(1,i);
        end 
        xLQR=xLQR(:,1:end-1);

        % Gráficas de los estados, salida y de la señal de control óptima
        % -----------Gráficas de las estimaciones del estado --------------
        figure; grid on; hold on
        plot(t,xPDIP(1,:)','-r','LineWidth', 3) 
        plot(t,xQuadProg(1,:)','-b','LineWidth', 2) 
        plot(t,xLQR(1,:),'--g','LineWidth', 1);
        legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
        title('Estado 1 - Velocidad angular N:'+string(N_HOR));
        box on; xlabel('Tiempo [s]'); ylabel('Velocidad Angular [rad/s]'); 
        grid on; box on;  hold off;

        figure; grid on; hold on;
        plot(t,xPDIP(2,:)','-r','LineWidth', 3) 
        plot(t,xQuadProg(2,:)','-b','LineWidth', 2) 
        plot(t,xLQR(2,:),'--g','LineWidth', 1);
        plot(t,r,':','LineWidth', 2);
        legend('MPC_{pdip}','MPC_{quadprog}','LQR', 'Referencia');  
        title('Estado 2 - Posición angular N:'+string(N_HOR)')
        ylim([-2,2]);
        box on; xlabel('Tiempo [s]'); ylabel(' Posisión Angular [rad]'); 
        grid on; box on;  hold off;
        
        % -----------Grafico de la entrada --------------------------------
        figure; grid on; hold on;
        plot(t,uPDIP,'-r','LineWidth', 3);
        plot(t,uQuadProg(1,:),'-b','LineWidth', 2);
        plot(t,uLQR,'--g','LineWidth', 1); 
        legend('MPC_{pdip}','MPC_{quadprog}','LQR');  
        title('Voltaje Entrada N:'+string(N_HOR)');
        ylim([-10,10]);
        ylabel('Voltaje [V]'); xlabel('Tiempo [s]'); 
        grid on; box on;  hold off;
    end
end


    




