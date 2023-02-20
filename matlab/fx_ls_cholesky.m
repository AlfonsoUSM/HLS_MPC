%% LINEAR SOLVER - CHOLESKY
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 09-02-2023
% Based on the work by Andrew Morrison https://github.com/morrisort/embeddedMPC/
% ===============================================================================

function [z] = myChol(B,v)
    n = size(B,1);
    for i=2:1:n
        for j=1:1:i-1
            B(i,j)=0;
        end
    end
    A=B;
    tol=n*eps;
    if A(1,1) <= tol
        A(1,1:n) = 0;
    else
        A(1,1:n) = A(1,1:n)/sqrt(A(1,1));
    end
    for j=2:n
      A(j,j:n) = A(j,j:n) - A(1:j-1,j)'*A(1:j-1,j:n);
      if A(j,j) <= tol
          A(j,j:n) = 0;
      else
          A(j,j:n) = A(j,j:n)/sqrt(A(j,j));
      end
    end
    suma=0; LT=A';
    y=zeros(n,1);
    for i=1:1:n
        for t=1:1:i-1
            suma=suma+LT(i,t)*y(t);
        end
        y(i,1)=(1/LT(i,i))*(v(i)-suma);
        suma=0;
    end
    suma=0;
    x=zeros(n,1);
    for i=0:1:n-1
        for t=1:1:i
            suma=suma+A(n-i,n-t+1)*x(n-t+1);
        end
        x(n-i,1)=(1/A(n-i,n-i))*(y(n-i)-suma);
        suma=0;
    end
    z =x;
end