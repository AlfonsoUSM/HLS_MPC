function [xinfy,uinfy] = fx_stationary(A,B,C,ref,Bpd,dis)
    n=size(A,1);
    m=size(B,2);
    T=[eye(n)-A,-B;C,zeros(m,m)];   
    if (nargin > 4)           
        bl=[Bpd*dis;ref];
    else            
        bl=[zeros(n,1);ref];
    end
    infy=T\bl; 
    xinfy=infy(1:n);         
    uinfy=infy(n+1:end);
end