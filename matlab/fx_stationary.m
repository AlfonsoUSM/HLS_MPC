function [xinfy,uinfy] = fx_stationary(A,B,C,ref,Bd,dis)
    n=size(A,1);
    m=size(B,2);
    S=[eye(n)-A,-B;C,zeros(m,m)];   
    if (nargin > 4)           
        bl=[Bd*dis;ref];
    else            
        bl=[zeros(n,1);ref];
    end
    infy=S\bl; 
    xinfy=infy(1:n);         
    uinfy=infy(n+1:end);
end