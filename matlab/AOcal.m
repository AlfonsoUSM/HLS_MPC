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
function [Acal,Ocal] = AOcal(A,B,N)
    Acal=[]; Ocal=[];  
    for i=1:1:N
        AB=[];
        Acal=[Acal;A^i];
        for j=1:1:i; AB=[AB,A^(i-j)*B]; end
        Fl=[AB,zeros(size(AB,1),N-i)];
        Ocal=[Ocal;Fl];
    end
end