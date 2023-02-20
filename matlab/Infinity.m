
% Cálculo "x" y "u" en estado estacionario
% Ad -> Matriz Ad de modelo espacio estado discreto
% Bd -> Matriz Bd de modelo espacio estado discreto
% Bpd -> Matriz Bpd de modelo espacio estado discreto
% Cd -> Matriz Cd de modelo espacio estado discreto
% per -> Perturbacion 
% yref -> Referencia de salida
function [xinfy, uinfy] = Infinity(Ad, Bd, Bpd, Cd, per, yref)
    H_k     = [ 0 0 1 0;  % Seleccionar vd y vq para hacer seguimiento
                0 0 0 1]; % de referencia
    % Asolve * infy =  Csolve
    Asolve  = [eye(length(Ad))-Ad          -Bd; 
                           H_k*Cd   zeros(2,2)];
    Csolve  = [[Bpd Bd]*per;
                       yref];
    
    % infy = Asolve \ Csolve
    infy  = linsolve(Asolve,Csolve);
    % infy = [xinfy; uinfy]
    xinfy = infy(1:size(Ad,1));         
    uinfy = infy(size(Ad,1)+1:end);
end
