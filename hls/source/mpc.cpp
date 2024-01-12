
#include "mpc.hpp"
#include "utils.hpp"

// ADMM Global Variables

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};


// Function Definitions

void mpc(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&u0)[M_SYS], int IT){
//#pragma HLS INTERFACE mode=ap_ctrl_none port=return
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=s_axilite port=x0
#pragma HLS INTERFACE mode=s_axilite port=u0
#pragma HLS INTERFACE mode=s_axilite port=IT
	// constraint c_hat = constraint(x0, r0)
#ifdef DENSE
	data_t g[M_QP];
	data_t q_hat[N_QP];
	data_t inf[N_SYS+M_SYS];
	mpc_dense_constraint(x0, r0, inf, q_hat, g);
#else
	data_t g[M_QP] = {0};
	data_t (&q_hat)[N_QP] = q;
	mpc_sparse_constraint(x0, g);
#endif
	// optimized theta = qp_solver(H,h_nau, C_hat, c_hat)
	qp_admm(q_hat, g, IT);
	for (int i=0; i<M_SYS ; i++){
#ifdef DENSE
		u0[i] = tk_admm[i] + inf[N_SYS+i];
#else
		u0[i] = tk_admm[(N_HOR*N_SYS+N_SYS + i)];
#endif
	}
	// estimate x1
	return;
}

#if defined DENSE
void mpc_dense_constraint(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&inf)[N_SYS+M_SYS], data_t (&q)[N_QP], data_t (&g)[M_QP]){
	// follow reference
	data_t ref[N_SYS+M_SYS];
	ref1: for (int i=0; i<N_SYS; i++){
		ref[i] = 0;
	}
	ref2: for (int i=0; i<P_SYS; i++){
		ref[i+N_SYS] = r0[i];
	}
	mvmult<(N_SYS+M_SYS),(N_SYS+P_SYS),data_t>(T_inv, ref, inf);
	data_t xnau[N_SYS];
	data_t xnau_max[N_SYS];
	data_t xnau_min[N_SYS];
	data_t unau_max[M_SYS];
	data_t unau_min[M_SYS];
	xinf: for (int i=0; i<N_SYS; i++){		// inf = [xinf; uinf]
		xnau[i] = x0[i] - inf[i];			// xnau = x - xinf
		xnau_max[i] = xmax[i] - inf[i];
		xnau_min[i] = xmin[i] - inf[i];
	}
	uinf: for (int i=0; i<M_SYS; i++){
		unau_max[i] = umax[i] - inf[N_SYS+i];
		unau_min[i] = umin[i] - inf[N_SYS+i];
	}
	// build vectors q and g
	vmmult<N_SYS,N_QP,data_t>(xnau, F, q);		// q = (xnau'*F)';
	data_t dnau_neg[N_QP], cnau[N_QP];
	ab1: for (int i=0; i<M_SYS; i++){
		ab2: for (int j=0; j<N_HOR; j++){
			cnau[j*M_SYS+i] = unau_max[i];
			dnau_neg[j*M_SYS+i] = -unau_min[i];
		}
	}
	data_t enau[N_SYS*N_HOR], fnau[N_SYS*N_HOR], temp[N_SYS*N_HOR];
	// f = [e-D*x0; D*x0-d];
	mvmult<(N_SYS*N_HOR),N_SYS,data_t>(D, xnau, temp);
	f1: for (int i=0; i<N_SYS; i++){
		f2: for (int j=0; j<N_HOR; j++){
			enau[j*N_SYS+i] = xnau_max[i];
			fnau[j*N_SYS+i] = xnau_min[i];
		}
	}
    // h = [f; b; -a];
	int i = 0;
	constraint1: for (int j=0; j<N_QP; j++){
		g[i] = cnau[j];
		i++;
	}
	constraint2: for (int j=0; j<N_QP; j++){
		g[i] = dnau_neg[j];
		i++;
	}
	constraint3: for (int j=0; j<(N_SYS*N_HOR); j++){
		g[i] = enau[j] - temp[j];
		i++;
	}
	constraint4: for (int j=0; j<(N_SYS*N_HOR); j++){
		g[i] = temp[j] - fnau[j];
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

void qp_admm(data_t (&q)[N_QP], data_t (&g)[M_QP], int IT){
	loop_admm: for(int i = 0; i < IT; i++){
#pragma HLS LOOP_TRIPCOUNT max=10
			data_t vx[M_QP];
			data_t temp[M_QP];
			admm_merge1:{	// vx = zk - g + uk;
				vsub<M_QP,data_t>(zk_admm, g, temp);
				vadd<M_QP,data_t>(temp, uk_admm, vx);
			}
			data_t temp1[N_QP], temp2[N_QP], temp3[N_QP];
			admm_merge2:{	// tk = R_inv * (-rho * G^T * vx - q);
				mvmult<N_QP,M_QP,data_t>(P, vx, temp1);
				vsub<N_QP,data_t>(temp1, q, temp2);
			}
			mvmult<N_QP,N_QP,data_t>(R_inv, temp2, tk_admm);
			data_t Gtk[M_QP], temp4[M_QP], temp5[M_QP], temp6[M_QP], temp7[M_QP];
			admm_merge:{	// zk = max{0, g - uk - G*tk};	uk = uk + (G*tk + zk - g);
				mvmult<M_QP,N_QP,data_t>(G, tk_admm, Gtk);		// Gtk = G*kt
				vsub<M_QP,data_t>(g, uk_admm, temp4);			// temp4 = g - uk
				vsub<M_QP,data_t>(temp4, Gtk, temp5);			// temp5 = (g - uk) - Gtk
				vadd<M_QP,data_t>(uk_admm, Gtk, temp6);			// temp6 = uk + Gtk
				max0<M_QP,data_t>(temp5, zk_admm);				// zk = max{0, g - uk - G*tk}
				vsub<M_QP,data_t>(temp6, g, temp7);				// temp7 = (uk + Gtk) - g
			}
				vadd<M_QP,data_t>(temp7, zk_admm, uk_admm);		// uk = ((uk + gtk) - g) + zk
			//}
	}
	return;
}
