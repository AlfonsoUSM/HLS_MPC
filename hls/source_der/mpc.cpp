
#include "mpc.hpp"
#include "utils.hpp"

// ADMM Global Variables

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};


// Function Definitions

void mpc(data_t (&x0)[N_SYS], data_t (&u0)[M_SYS], int IT){
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=s_axilite port=x0
#pragma HLS INTERFACE mode=s_axilite port=u0
#pragma HLS INTERFACE mode=s_axilite port=IT
	// constraint c_hat = constraint(x0, r0)
#if defined DENSE
	data_t h[M_QP];
	data_t q_hat[N_QP];
	mpc_dense_constraint(x0, q_hat, h);
#else
	data_t h[M_QP] = {0};
	data_t (&q_hat)[N_QP] = q;
	mpc_sparse_constraint(x0, h);
#endif
	// optimized theta = qp_solver(H,h_nau, C_hat, c_hat)
	qp_admm(q_hat, h, IT);
	for (int i=0; i<M_SYS ; i++){
#ifdef DENSE
		u0[i] = tk_admm[i];
#else
		u0[i] = tk_admm[(N_HOR*N_SYS+N_SYS + i)];
#endif
	}
	// estimate x1
	return;
}

#if defined DENSE
void mpc_dense_constraint(data_t (&x0)[N_SYS], data_t (&q)[N_QP], data_t (&h)[M_QP]){
	// follow reference currently not implemented
	vmmult<N_SYS,N_QP,data_t>(x0, G, q);		// q = (x0'*G)';
	data_t temp[5*N_HOR], f1[5*N_HOR], f2[5*N_HOR];
	mvmult<(5*N_HOR),N_SYS,data_t>(D, x0, temp);
	vsub<(5*N_HOR),data_t>(e, temp, f1);
	vsub<(5*N_HOR),data_t>(temp, d, f2);	// f = [e-D*x0; D*x0-d];
    // h = [f; b; -a];
	int i = 0;
	constraint1: for (int j=0; j<(5*N_HOR); j++){
		h[i] = f1[j];
		i++;
	}
	constraint2: for (int j=0; j<(5*N_HOR); j++){
		h[i] = f2[j];
		i++;
	}
	constraint3: for (int j=0; j<(5*N_HOR); j++){
		h[i] = b[j];
		i++;
	}
	constraint4: for (int j=0; j<(5*N_HOR); j++){
		h[i] = a_neg[j];
		i++;
	}
	return;
}

#else
void mpc_sparse_constraint(data_t (&x0)[N_SYS], data_t (&h)[M_QP]){
	// follow reference currently not implemented
	// h = [g; f; -f], where f = [-x0; 0; 0; 0...] is (N_HOR+1)*N_SYS
	constraint1: for (int i=0; i<(2*N_QP); i++){
		h[i] = g[i];
	}
	constraint2: for (int i=0; i<N_SYS; i++){
		h[(2*N_QP + i)] = -x0[i];
	}
	constraint3: for (int i=0; i<N_SYS; i++){
		h[(2*N_QP + N_HOR*N_SYS + N_SYS + i)] = x0[i];
	}
	return;
}
#endif

void qp_admm(data_t (&q)[N_QP], data_t (&h)[M_QP], int IT){
	loop_admm: for(int i = 0; i < IT; i++){
#pragma HLS LOOP_TRIPCOUNT max=10
			data_t vx[M_QP];
			data_t temp[M_QP];
			admm_merge1:{
	// vx = zk - h + uk;
				vsub<M_QP,data_t>(zk_admm, h, temp);
				vadd<M_QP,data_t>(temp, uk_admm, vx);
			}
			data_t temp1[N_QP], temp2[N_QP], temp3[N_QP];
			admm_merge2:{
	// tk = R_inv * (-rho * H^T * vx - q);
				mvmult<N_QP,M_QP,data_t>(W, vx, temp1);
				vsub<N_QP,data_t>(temp1, q, temp2);
			}
			mvmult<N_QP,N_QP,data_t>(R_inv, temp2, tk_admm);
			data_t Htk[M_QP], temp4[M_QP], temp5[M_QP], temp6[M_QP], temp7[M_QP];
			admm_merge:{
	// zk = max{0, h - uk - H*tk};	uk = uk + (H*tk + zk - h);
				mvmult<M_QP,N_QP,data_t>(H, tk_admm, Htk);		// Htk = H*kt
				vsub<M_QP,data_t>(h, uk_admm, temp4);			// temp4 = h - uk
				vsub<M_QP,data_t>(temp4, Htk, temp5);			// temp5 = (h - uk) - Htk
				vadd<M_QP,data_t>(uk_admm, Htk, temp6);			// temp6 = uk + Htk
				max0<M_QP,data_t>(temp5, zk_admm);				// zk = max{0, h - uk - H*tk}
				vsub<M_QP,data_t>(temp6, h, temp7);				// temp7 = (uk + Htk) - h
			}
				vadd<M_QP,data_t>(temp7, zk_admm, uk_admm);		// uk = ((uk + Htk) - h) + zk
			//}
	}
	return;
}
