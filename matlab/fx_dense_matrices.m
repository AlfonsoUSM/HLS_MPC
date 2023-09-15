% Función para generar las matrices de Acal*x_0 + Ocal*vec{u}
% Acal = [   A
%            A^2
%            .
%            .
%            A^(N)]
%        
% Ocal = [   B
%            AB          B
%            .           .        .
%            A^(N-1)B    A^(N-2)B . . . B]       
function [D,E] = fx_dense_matrices(A,B,N)
    D=[]; E=[];  
    for i=1:1:N
        AB=[];
        D=[D;A^i];
        for j=1:1:i 
            AB=[AB,A^(i-j)*B]; 
        end
        Fl=[AB,zeros(size(AB,1),N-i)];
        E=[E;Fl];
    end
end