%% Poligono circunscrito en una circunferencia
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 17-01-2024
% Based on the work by Reinier López
% ===============================================================================


function [H, max] = fx_poliedro(R, N)
% Entrega poligono con N lados circunscrito en la circunferencia de radio R.
% Para cada recta del conjunto H se tiene H(i,:)*[x; y] ~ R*max[i]
% N must be an even number
    
    % Inicializar matrices para almacenar los coeficientes A y B, y el vector C
    A = zeros(N, 2);
    max_aux = zeros(N, 1);
    % Inicializar matrices para almacenar las coordenadas
    coordenadas_x = zeros(1, N);
    coordenadas_y = zeros(1, N);
    
    % Calcular coordenadas de los puntos
    angle = 2*pi/N;         % Calculo del angulo
    for i = 1:N
        theta = i * angle;
        coordenadas_x(i) = R * cos(theta);
        coordenadas_y(i) = R * sin(theta);
    end
    
    % Calcular coeficientes H, [x, y] y max para cada ecuación de la recta
    for i = 1 : (N)
        x1 = coordenadas_x(i);
        y1 = coordenadas_y(i);
        if i < N
            x2 = coordenadas_x(i + 1);  % Circular: conecta el último con el primero
            y2 = coordenadas_y(i + 1);
        end
        if i == N
            x2 = coordenadas_x(1);  % Circular: conecta el último con el primero
            y2 = coordenadas_y(1);
        end
        % Calcular la pendiente y el término independiente
        m = (y2 - y1) / (x2 - x1);
        n = y1 - m * x1;
        A(i, :) = [-m, 1];          % Asignar los coeficientes A y B
        max_aux(i) = -n/R;          % Asignar el vector max
    end

    max = zeros(N/2,1);
    H = zeros(N/2,2);
    j = 1;
    for i = 1:N
        if max_aux(i) > 0
            max(j) = max_aux(i);
            H(j,:) = A(i, :);
            j = j+1;
        end
    end
end