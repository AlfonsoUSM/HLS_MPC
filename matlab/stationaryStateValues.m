function [xinfy,uinfy] = stationaryStateValues(A,B,C,yref)
    n=size(A,1);
    m=size(B,2);
    p=size(C,1);
    L=[A-eye(n),B;C,zeros(p,m)];              
    bl=[zeros(n,1);yref];
    infy=L\bl; 
    xinfy=infy(1:n);         
    uinfy=infy(n+1:end);
end