
#include "mpc.hpp"
#include "utils.hpp"

// ADMM Global Variables

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};


// Function Definitions

void mpc(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&d0)[D_SYS], data_t (&u0)[M_SYS], int IT){
#pragma HLS INTERFACE mode=s_axilite port=return
#pragma HLS INTERFACE mode=s_axilite port=x0
#pragma HLS INTERFACE mode=s_axilite port=u0
#pragma HLS INTERFACE mode=s_axilite port=IT
	// constraint g_nau, q_nau = constraint(x0, r0, d0)
#if defined DENSE
	data_t g_nau[M_QP];
	data_t q_nau[N_QP];
	data_t inf[N_SYS+M_SYS];
	mpc_dense_constraint(x0, r0, d0, inf, q_nau, g_nau);
#else
	data_t h[M_QP] = {0};
	data_t (&q_hatnau)[N_QP] = q;
	mpc_sparse_constraint(x0, h);
#endif
	// optimized theta = qp_solver(Q,q_nau, G_nau, g_nau)
	qp_admm(q_nau, g_nau, IT);
	for (int i=0; i<M_SYS ; i++){
#ifdef DENSE
		u0[i] = tk_admm[i] + inf[N_SYS+i];
#else
		u0[i] = tk_admm[(N_HOR*N_SYS+N_SYS + i)];
#endif
	}
	return;
}

#if defined DENSE
void mpc_dense_constraint(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&d0)[D_SYS], data_t (&inf)[N_SYS+M_SYS], data_t (&q_nau)[N_QP], data_t (&g_nau)[M_QP]){
	// follow reference currently not implemented
	data_t Bpdd0[N_SYS];
	mvmult<N_SYS,D_SYS,data_t>(Bpd, d0, Bpdd0);
	data_t ref[N_SYS+M_SYS];
	ref1: for (int i=0; i<N_SYS; i++){
		ref[i] = Bpdd0[i];
	}
	ref2: for (int i=0; i<P_SYS; i++){
		ref[i+N_SYS] = r0[i];
	}
	mvmult<(N_SYS+M_SYS),(N_SYS+M_SYS),data_t>(T_inv, ref, inf);
	data_t uinf[M_SYS];
	data_t xinf[N_SYS];
	data_t xnau[N_SYS];
	xinf: for (int i=0; i<N_SYS; i++){
		xinf[i] = inf[i]; 					// inf = [xinf; uinf]
		xnau[i] = x0[i] - inf[i];			// xnau = x0 - xinf
	}
	uinf: for (int i=0; i<M_SYS; i++){
		uinf[i] = inf[N_SYS+i];
	}
	// build vectors q and g
	vmmult<N_SYS,N_QP,data_t>(xnau, F, q_nau);		// q = (xnau'*F)';
	data_t Huinf[A_SYS];
	data_t Ixinf[B_SYS];
	mvmult<A_SYS,M_SYS,data_t>(H, uinf, Huinf);
	mvmult<B_SYS,N_SYS,data_t>(I, xinf, Ixinf);
	data_t WDx[B_SYS*N_HOR];
	mvmult<B_SYS*N_HOR,N_SYS,data_t>(WD, xnau, WDx);
	data_t unau_max[A_SYS];
	data_t unau_min[A_SYS];
	data_t xnau_max[B_SYS];
	data_t xnau_min[B_SYS];
	vsub<A_SYS, data_t>(umax, Huinf, unau_max);
	vsub<A_SYS, data_t>(Huinf, umin, unau_min);
	vsub<B_SYS, data_t>(xmax, Ixinf, xnau_max);
	vsub<B_SYS, data_t>(Ixinf, xmin, xnau_min);
	// g = [c; d; e-WDx; f+WDx];
	constraint1: for (int i=0; i<A_SYS; i++){
		for (int j=0; j<N_HOR; j++){
			g_nau[A_SYS*j+i] = unau_max[i];
		}
		for (int j=N_HOR; j<(2*N_HOR); j++){
			g_nau[A_SYS*j+i] = unau_min[i];
		}
	}
	int k = 2*N_HOR*A_SYS;
	constraint2: for (int i=0; i<B_SYS; i++){
		for (int j=0; j<N_HOR; j++){
			g_nau[k+B_SYS*j+i] = xnau_max[i]-WDx[B_SYS*j+i];
		}
		for (int j=0; j<N_HOR; j++){
			g_nau[k+B_SYS*(j+N_HOR)+i] = xnau_min[i]+WDx[B_SYS*j+i];
		}
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
			admm_merge:{	// zk = max{0, g - uk - H*tk};	uk = uk + (G*tk + zk - g);
				mvmult<M_QP,N_QP,data_t>(G, tk_admm, Gtk);		// Gtk = G*kt
				vsub<M_QP,data_t>(g, uk_admm, temp4);			// temp4 = g - uk
				vsub<M_QP,data_t>(temp4, Gtk, temp5);			// temp5 = (g - uk) - Gtk
				vadd<M_QP,data_t>(uk_admm, Gtk, temp6);			// temp6 = uk + Gtk
				max0<M_QP,data_t>(temp5, zk_admm);				// zk = max{0, g - uk - G*tk}
				vsub<M_QP,data_t>(temp6, g, temp7);				// temp7 = (uk + Gtk) - g
			}
				vadd<M_QP,data_t>(temp7, zk_admm, uk_admm);		// uk = ((uk + Gtk) - g) + zk
			//}
	}
	return;
}
