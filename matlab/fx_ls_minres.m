%% LINEAR SOLVER - MINRES
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 09-02-2023
% Based on the work by Andrew Morrison https://github.com/morrisort/embeddedMPC/
% ===============================================================================

function [x,iter]=fx_ls_minres(A,b,x0,n_max,tol)
% computes the solution of Ax=b using the MINRES scheme. (book)
% xO: start vector 
% n_max: maximal number of iterations
% tol: desired accuracy of the residual
% initialize
    i=1; N=length(x0);                              
    v=zeros(N,1);                                    
    v_hat=b-A*x0;                          
    beta=norm(v_hat);                               
    c=1;                                            
    c_old=1;                                         
    s_old=0;                                         
    s=0;                                             
    w=zeros(N,1);                                    
    w_old=w;                                         
    eta=beta;                                        
    x=x0;                                           
    norm_rMR=beta;                                   
    norm_r0=beta;
    iter=i;
    while (i < n_max+1) && (norm_rMR/norm_r0 > tol)
        i=i+1;  
        iter=i;
        %%Lanczos
        v_old=v;                                            
        v=v_hat/beta;                                       
        Av=A*v;                                     
        alpha=v'*Av;                                        
        v_hat=Av-alpha*v-beta*v_old;                         
        beta_old=beta;                                       
        beta=norm(v_hat);   % norm(v)=sqrt(v'*v)
        %%QR factorization
        c_oold=c_old;                                       
        c_old=c;                                            
        s_oold=s_old;                                        
        s_old=s;                                            
        r1_hat=c_old*alpha-c_oold*s_old*beta_old;           
        r1 =sqrt(r1_hat^2+beta^2);                          
        r2 =s_old*alpha+c_oold*c_old*beta_old;              
        r3 =s_oold*beta_old;                                
        %%Givens rotation
        c=r1_hat/r1;                                        
        s=beta/r1;                                          
        %%update
        w_oold=w_old;                                       
        w_old=w;                                           
        w=(v-r3*w_oold-r2*w_old)/r1;                        
        x=x+c*eta*w;                                    
        norm_rMR=norm_rMR*abs(s);                           
        eta=-s*eta;                                         
    end 
%     if abs(c) > eps
%         xOR =x-s*eta*w/c;                                 
%         norm_rOR=norm_rMR/abs(c);                           
%     end
    return;
end