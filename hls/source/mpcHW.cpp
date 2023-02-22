#include"mpcHW.h"

void mpcHW(elem x_axi[N_SYS], elem r_axi[P_SYS], elem xmin_axi[N_SYS], elem xmax_axi[N_SYS], elem umin_axi[P_SYS], elem umax_axi[P_SYS], elem Acal_axi[N_SYS*N_HOR][N_SYS], elem AcalQOcal_axi[N_SYS][N_HOR], elem H_axi[N_HOR][N_HOR], elem M_hat_axi[I_AUX][N_HOR], elem L_invLast_axi[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u_star_axi[M_SYS]){
#pragma HLS INTERFACE s_axilite port=iterPDIP
#pragma HLS INTERFACE s_axilite port=iterLS

#pragma HLS INTERFACE m_axi port=x_axi
#pragma HLS INTERFACE m_axi port=r_axi
#pragma HLS INTERFACE m_axi port=xmin_axi
#pragma HLS INTERFACE m_axi port=xmax_axi
#pragma HLS INTERFACE m_axi port=umin_axi
#pragma HLS INTERFACE m_axi port=umax_axi
#pragma HLS INTERFACE m_axi port=Acal_axi
#pragma HLS INTERFACE m_axi port=AcalQOcal_axi
#pragma HLS INTERFACE m_axi port=H_axi
#pragma HLS INTERFACE m_axi port=M_hat_axi
#pragma HLS INTERFACE m_axi port=L_invLast_axi
#pragma HLS INTERFACE m_axi port=u_star_axi



	elem x[N_SYS];
	elem r[P_SYS];
	elem xmin[N_SYS];
	elem xmax[N_SYS];
	elem umin[P_SYS];
	elem umax[P_SYS];
	elem Acal[N_SYS*N_HOR][N_SYS];
	elem AcalQOcal[N_SYS][N_HOR];
	elem H[N_HOR][N_HOR];
	elem M_hat[I_AUX][N_HOR];
	elem L_invLast[N_SYS+P_SYS];
	elem u_star[P_SYS];

	// store values in local memory
	for(int i=0; i<N_SYS ; i++){
		x[i] = x_axi[i];
	}
	for(int i=0; i<P_SYS ; i++){
		r[i] = r_axi[i];
	}
	for(int i=0; i<N_SYS ; i++){
		xmin[i] = xmin_axi[i];
	}
	for(int i=0; i<N_SYS ; i++){
		xmax[i] = xmax_axi[i];
	}
	for(int i=0; i<P_SYS ; i++){
		umin[i] = umin_axi[i];
	}
	for(int i=0; i<P_SYS ; i++){
		umax[i] = umax_axi[i];
	}
	for(int r=0; r< (N_SYS*N_HOR); r++){
		for(int c=0; c< N_SYS; c++){
			Acal[r][c] = Acal_axi[r][c];
		}
	}
	for(int r=0; r< N_SYS; r++){
		for(int c=0; c< N_HOR; c++){
			AcalQOcal[r][c] = AcalQOcal_axi[r][c];
		}
	}
	for(int r=0; r< N_HOR; r++){
		for(int c=0; c< N_HOR; c++){
			H[r][c] = H_axi[r][c];
		}
	}
	for(int r=0; r< I_AUX; r++){
		for(int c=0; c< N_HOR; c++){
			M_hat[r][c] = M_hat_axi[r][c];
		}
	}
	for(int i=0; i< (N_SYS+P_SYS); i++){
		L_invLast[i] = L_invLast_axi[i];
	}

	elem c_hat[I_AUX];
	elem h_tilde[N_HOR];
#pragma HLS ARRAY_PARTITION variable=c_hat complete
#pragma HLS ARRAY_PARTITION variable=x complete
#pragma HLS ARRAY_PARTITION variable=r complete
#pragma HLS ARRAY_PARTITION variable=L_invLast complete
	elem x_inf[N_SYS];
	elem u_inf[M_SYS];

	// MATLAB: [x_inf,u_inf] = stationaryStateValues(A,B,C,r);
	set_x_inf:
	for (int i=0; i<N_SYS; i++){
#pragma HLS UNROLL
		x_inf[i] = L_invLast[i]*r[0];
	}
	set_u_inf:
	for (int i=0; i<M_SYS; i++){
#pragma HLS UNROLL
		u_inf[i] = L_invLast[i+N_SYS]*r[0];
	}

	elem a_tilde[N_HOR];
	elem b_tilde[N_HOR];
#pragma HLS ARRAY_PARTITION variable=a_tilde complete
#pragma HLS ARRAY_PARTITION variable=b_tilde complete
	// MATLAB: a_tilde=ones(N,1)*(umin-u_inf); 
	// MATLAB: b_tilde=ones(N,1)*(umax-u_inf);
	set_a:
	for (int i=0; i<N_HOR; i++){
#pragma HLS UNROLL
		a_tilde[i] = umin[0] - u_inf[0];
	}
	set_b:
	for (int i=0; i<N_HOR; i++){
#pragma HLS UNROLL
		b_tilde[i] = umax[0] - u_inf[0];
	}
	elem x_tilde[N_SYS];
	// MATLAB: x_tilde=x-x_inf;
	set_x_tilde:
	for (int i=0; i<N_SYS; i++){
#pragma HLS UNROLL
		x_tilde[i] = x[i] - x_inf[i];
	}

	// MATLAB: h_tilde=(2*x_tilde'*Acal'*Q*Ocal)';
	get_h_tilde(AcalQOcal, x_tilde, h_tilde);


	// MATLAB: c_tilde=[(Xmax-Acal*x_tilde);-(Xmin-Acal*x_tilde)];
	// MATLAB: c_hat=[c_tilde;b_tilde;-a_tilde];
	get_c_hat(Acal, xmax, xmin, x_inf, x_tilde, a_tilde, b_tilde, c_hat);
	elem tol = 1e-9;
	elem u_tilde_star[N_HOR];

	// MATLAB: [u_tilde_star,~,~,~,saveMat]=pdip(H,h_tilde,M_hat,c_hat,iterPDIP,iterMINRES,linSol,tol,saveMat);
	pdipHW(H, h_tilde, M_hat, c_hat, u_tilde_star, iterPDIP, iterLS, tol);

	// MATLAB: u_star=u_tilde_star(1)+u_inf;
	for (int i=0; i<P_SYS; i++){
		u_star[i]=u_tilde_star[i]+u_inf[i];
	}
	// AXI output
	for (int i=0; i<P_SYS; i++){
		u_star_axi[i]=u_star[i];
	}
}

void get_h_tilde(elem AcalQOcal[N_SYS][N_HOR], elem x_tilde[N_SYS], elem h_tilde[N_HOR]){
#pragma HLS INLINE off
	set_h_tilde:
	for (int i=0; i<N_HOR; i++){
#pragma HLS pipeline
		h_tilde[i] = 0;
		for (int j=0; j<N_SYS; j++){
			h_tilde[i] += 2*x_tilde[j]*AcalQOcal[j][i];
		}
	}
}
void get_c_hat(elem Acal[N_SYS*N_HOR][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_inf[N_SYS], elem x_tilde[N_SYS], elem a_tilde[N_HOR], elem b_tilde[N_HOR], elem c_hat[6*N_HOR]){
	set_c_hat_c1:
	for (int i=0; i<2*N_HOR; i++){
#pragma HLS pipeline
		c_hat[i] =   xmax[i%N_SYS]- x_inf[i%N_SYS];
		for (int c=0; c<N_SYS; c++){
			c_hat[i] -= Acal[i%(N_SYS*N_HOR)][c]*x_tilde[c];
		}
	}
	set_c_hat_c2:
	for (int i=0; i<2*N_HOR; i++){
#pragma HLS pipeline
		c_hat[i+2*N_HOR] =   - (xmin[i%N_SYS]- x_inf[i%N_SYS]);
		for (int c=0; c<N_SYS; c++){
			c_hat[i+2*N_HOR] += Acal[i%(N_SYS*N_HOR)][c]*x_tilde[c];
		}
	}
	set_cx_b:
	for (int i=0; i<N_HOR; i++){
#pragma HLS pipeline
		c_hat[i+4*N_HOR] = b_tilde[i];
	}
	set_cx_a:
	for (int i=0; i<N_HOR; i++){
#pragma HLS pipeline
		c_hat[i+5*N_HOR] = -a_tilde[i];
	}
}