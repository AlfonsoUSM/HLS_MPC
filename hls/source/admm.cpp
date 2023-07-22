
#include "admm.hpp"
#include "utils.hpp"

data_t Rho = 2;
data_t R_inv[N_QP][N_QP] = {2};
data_t RhoMt_neg[N_QP][M_QP] = {2};

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};

void qp_admm(data_t (&c_qp)[M_QP]){
	loop_admm: for(int i = 0; i < N_IT; i++){
// #pragma HLS UNROLL factor=2
	        tk_update(c_qp, tk_admm, zk_admm, uk_admm);
	        //zk_update(c_qp, tk_admm, zk_admm, uk_admm);
	        //uk_update(c_qp, tk_admm, zk_admm, uk_admm);
	        zk_uk_update(c_qp, tk_admm, zk_admm, uk_admm);
	    }
	return;
}


void tk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t vx[M_QP];
	data_t temp[M_QP];
	// vx = zk - c + uk;
	tk_update_merge1:{
		vsub<M_QP,data_t>(zk, c_qp, temp);
		vadd<M_QP,data_t>(temp, uk, vx);
	}
    // tk = R_inv * (-rho * M^T * vx - h);
    data_t temp1[N_QP], temp2[N_QP], temp3[N_QP];
    tk_update_merge2:{
    	mvmult<N_QP,M_QP,data_t>(RhoMt_neg, vx, temp1);
    	vsub<N_QP,data_t>(temp1, h_qp, temp2);
    }
    mvmult<N_QP,N_QP,data_t>(R_inv, temp2, tk);
    return;
}

void zk_uk_update (data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
//#pragma HLS ARRAY_PARTITION variable=tk type=complete
	data_t Ctk[M_QP], temp[M_QP], temp1[M_QP], temp2[M_QP], temp3[M_QP];
	zk_uk_update_merge:{
    //  zk = max{0, c - uk - C*tk};
    // uk = uk + (C*tk + zk - c);
	mvmult<M_QP,N_QP,data_t>(C_qp, tk, Ctk);	// Ctk = C*kt
	vsub<M_QP,data_t>(c_qp, uk, temp);			// temp = c - uk
	vsub<M_QP,data_t>(temp, Ctk, temp1);		// temp1 = (c - uk) - Ctk
	vadd<M_QP,data_t>(uk, Ctk, temp2);			// temp2 = uk + Ctk
	max0<M_QP,data_t>(temp1, zk);				// zk
	vsub<M_QP,data_t>(temp2, c_qp, temp3);		// temp3 = (uk + Ctk) - c
	vadd<M_QP,data_t>(temp3, zk, uk);			// uk = ((uk + Ctk) - c) + zk
	}
	return;
}

void zk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t temp[M_QP], temp1[M_QP], temp2[M_QP];
    {
//#pragma HLS LOOP_MERGE
    //  zk = max{0, c - uk - M*tk};
	vsub<M_QP,data_t>(c_qp, uk, temp);
	mvmult<M_QP,N_QP,data_t>(C_qp, tk, temp1);
	vsub<M_QP,data_t>(temp, temp1, temp2);
	max0<M_QP,data_t>(temp2, zk);
    }
	return;
}

void uk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t temp[M_QP], temp1[M_QP], temp2[M_QP];
    {
//#pragma HLS LOOP_MERGE
    // uk = uk + (M*tk + zk - c);
	vsub<M_QP,data_t>(zk, c_qp, temp);
	vadd<M_QP,data_t>(uk, temp, temp1);
	mvmult<M_QP,N_QP,data_t>(C_qp, tk, temp2);
	vadd<M_QP,data_t>(temp1, temp2, uk);
    }
	return;
}

