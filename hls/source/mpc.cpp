
#include "mpc.hpp"
#include "utils.hpp"

// ADMM Global Variables

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};


// Function Definitions

void mpc(data_t (&x0)[N_SYS], data_t (&u0)[M_SYS]){
	data_t h[M_QP] = {0};
	data_t r0[P_SYS] = {0};
	// constraint c_hat = constraint(x0, r0)
#if defined DENSE
	mpc_dense_constraint(r0, x0, h);
#else
	mpc_sparse_constraint(r0, x0, h);
#endif
	// optimized theta = qp_solver(H,h_nau, C_hat, c_hat)
	qp_admm(h);
	for (int i=0; i<M_SYS ; i++){
		u0[i] = tk_admm[(N_HOR*N_SYS+N_SYS + i)];
	}
	// estimate x1
	return;
}

void mpc_dense_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&h)[M_QP]){
	// follow reference currently not implemented
	return;
}

void mpc_sparse_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&h)[M_QP]){
	// follow reference currently not implemented
	// c_hat = [g; f; -f], where f = [-x0; 0; 0; 0...] is (N_HOR+1)*N_SYS
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

void qp_admm(data_t (&h)[M_QP]){
	loop_admm: for(int i = 0; i < N_IT; i++){
			data_t vx[M_QP];
			data_t temp[M_QP];
			// vx = zk - h + uk;
			admm_merge1:{
				vsub<M_QP,data_t>(zk_admm, h, temp);
				vadd<M_QP,data_t>(temp, uk_admm, vx);
			}
			// tk = R_inv * (-rho * H^T * vx - q);
			data_t temp1[N_QP], temp2[N_QP], temp3[N_QP];
			admm_merge2:{
				mvmult<N_QP,M_QP,data_t>(RhoHt_neg, vx, temp1);
				vsub<N_QP,data_t>(temp1, q, temp2);
			}
			mvmult<N_QP,N_QP,data_t>(R_inv, temp2, tk_admm);
		    //  zk = max{0, h - uk - H*tk};
		    // uk = uk + (H*tk + zk - h);
			data_t Htk[M_QP], temp4[M_QP], temp5[M_QP], temp6[M_QP], temp7[M_QP];
			admm_merge:{
				mvmult<M_QP,N_QP,data_t>(H, tk_admm, Htk);	// Htk = H*kt
				vsub<M_QP,data_t>(h, uk_admm, temp4);			// temp = h - uk
				vsub<M_QP,data_t>(temp4, Htk, temp5);		// temp1 = (h - uk) - Htk
				vadd<M_QP,data_t>(uk_admm, Htk, temp6);			// temp2 = uk + Htk
				max0<M_QP,data_t>(temp5, zk_admm);				// zk
				vsub<M_QP,data_t>(temp6, h, temp7);		// temp3 = (uk + Htk) - c
				vadd<M_QP,data_t>(temp7, zk_admm, uk_admm);			// uk = ((uk + Htk) - h) + zk
			}

	    }
	return;
}
