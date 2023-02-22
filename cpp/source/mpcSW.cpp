#include"mpcSW.h"

void mpcSW(elem x[N_SYS], elem yref[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[P_SYS], elem umax[P_SYS], elem Acal[N_SYS*N_HOR][N_SYS], elem AcalQOcal[N_SYS][N_HOR], elem H[N_HOR][N_HOR], elem M_hat[I_AUX][N_HOR], elem L_invLast[N_SYS+P_SYS], int iterPDIP, int iterLS, elem tol, elem u_star[M_SYS]){

	elem c_hat[6*N_HOR];
    elem h_tilde[N_HOR];
	elem x_inf[N_SYS];
	elem u_inf[P_SYS];
	elem a_tilde[N_HOR];
	elem b_tilde[N_HOR];
	elem x_tilde[N_SYS];
	elem u_tilde[N_HOR];

	// MATLAB: [x_inf,u_inf] = stationaryStateValues(A,B,C,r);
	for (int i=0; i<N_SYS; i++){
		x_inf[i] = L_invLast[i]*yref[0];
	}
	for (int i=0; i<P_SYS; i++){
		u_inf[i] = L_invLast[i+N_SYS]*yref[0];
	}
	
	// MATLAB: a_tilde=ones(N,1)*(umin-u_inf); 
	for (int i=0; i<N_HOR; i++){
		a_tilde[i] = umin[0] - u_inf[0];
	}
	// MATLAB: b_tilde=ones(N,1)*(umax-u_inf);	
	for (int i=0; i<N_HOR; i++){
		b_tilde[i] = umax[0] - u_inf[0];
	}
	
	// MATLAB: x_tilde=x-x_inf;
	for (int i=0; i<N_SYS; i++){
		x_tilde[i] = x[i] - x_inf[i];
	}
	
	// MATLAB: h_tilde=(2*x_tilde'*Acal'*Q*Ocal)';
	// AclaQOcal =Acal'*Q*Ocal
    get_h_tilde(AcalQOcal, x_tilde, h_tilde);

	// MATLAB: c_tilde=[(Xmax-Acal*x_tilde);-(Xmin-Acal*x_tilde)];
	// MATLAB: c_hat=[c_tilde;b_tilde;-a_tilde];
	get_c_hat(Acal, xmax, xmin, x_inf, x_tilde, a_tilde, b_tilde, c_hat);

	// MATLAB: [u_tilde_star,~,~,~,saveMat]=pdip(H,h_tilde,M_hat,c_hat,iterPDIP,iterMINRES,linSol,tol,saveMat);
    pdipSW(H, h_tilde, M_hat, c_hat, u_tilde, iterPDIP, iterLS, tol);

	// MATLAB: u_star=u_tilde_star(1)+u_inf;
	u_star[0]=u_tilde[0]+u_inf[0];
}

void get_h_tilde(elem AcalQOcal[N_SYS][N_HOR], elem x_tilde[N_SYS], elem h[N_HOR]){
	for (int i=0; i<N_HOR; i++){
        h[i] = 0;
		for (int j=0; j<N_SYS; j++){
            h[i] += 2*x_tilde[j]*AcalQOcal[j][i];
		}
	}
}
void get_c_hat(elem Acal[N_SYS*N_HOR][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_inf[N_SYS], elem x_tilde[N_SYS], elem a[N_HOR], elem b[N_HOR], elem c_hat[I_AUX]){

	for (int i=0; i<2*N_HOR; i++){
		c_hat[i] =   xmax[i%N_SYS]- x_inf[i%N_SYS];
		for (int c=0; c<N_SYS; c++){
			c_hat[i] -= Acal[i%(N_SYS*N_HOR)][c]*x_tilde[c];
		}
	}
	
	for (int i=0; i<2*N_HOR; i++){
		c_hat[i+2*N_HOR] =   - (xmin[i%N_SYS]- x_inf[i%N_SYS]);
		for (int c=0; c<N_SYS; c++){
			c_hat[i+2*N_HOR] += Acal[i%(N_SYS*N_HOR)][c]*x_tilde[c];
		}
	}
	
	for (int i=0; i<N_HOR; i++){
		c_hat[i+4*N_HOR] = b[i];
	}
	
	for (int i=0; i<N_HOR; i++){
		c_hat[i+5*N_HOR] = -a[i];
	}
}
