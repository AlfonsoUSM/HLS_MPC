function [xinfy,uinfy] = fx_stationary(A,B,C,ref)
    n=size(A,1);
    m=size(B,2);
    S=[eye(n)-A,-B;C,zeros(m,m)];              
    bl=[zeros(n,1);ref];
    infy=S\bl; 
    xinfy=infy(1:n);         
    uinfy=infy(n+1:end);
end