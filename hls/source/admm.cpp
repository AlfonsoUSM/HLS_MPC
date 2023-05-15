
#include "admm.hpp"
#include "utils.hpp"

void qp_admm(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
#pragma HLS ARRAY_PARTITION variable=tk type=complete
	loop_admm: for(int i = 0; i < N_IT; i++){
// #pragma HLS UNROLL factor=2
	        tk_update(c_qp, tk, zk, uk);
	        //zk_update(c_qp, tk, zk, uk);
	        //uk_update(c_qp, tk, zk, uk);
	        zk_uk_update(c_qp, tk, zk, uk);
	    }
	return;
}


void tk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t vx[M_QP];
	data_t temp[M_QP];
	// vx = zk - c + uk;
	{
#pragma HLS LOOP_MERGE force
		vsub<M_QP,data_t>(zk, c_qp, temp);
		vadd<M_QP,data_t>(temp, uk, vx);
	}
    // tk = R_inv * (-rho * A^T * vx - h);
    data_t temp1[N_QP], temp2[N_QP], temp3[N_QP];
    {
#pragma HLS LOOP_MERGE force
    	mvmult<N_QP,M_QP,data_t>(RhoMt_neg, vx, temp1);
    	vsub<N_QP,data_t>(temp1, h_qp, temp2);
    }
    mvmult<N_QP,N_QP,data_t>(R_inv, temp2, tk);
    return;
}

void zk_uk_update (data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
//#pragma HLS ARRAY_PARTITION variable=tk type=complete
	data_t Mtk[M_QP], temp[M_QP], temp1[M_QP], temp2[M_QP], temp3[M_QP];
	{
#pragma HLS LOOP_MERGE
    //  zk = max{0, c - uk - M*tk};
    // uk = uk + (M*tk + zk - c);
	vsub<M_QP,data_t>(c_qp, uk, temp);			// temp = c - uk
	vsub<M_QP,data_t>(zk, c_qp, temp1);			// temp1 = zk - c
	vadd<M_QP,data_t>(uk, temp1, temp2);		// temp2 = uk + (zk - c)
	mvmult<M_QP,N_QP,data_t>(M_qp, tk, Mtk);	// M*kt
	vsub<M_QP,data_t>(temp, Mtk, temp3);		// temp2 = (c - uk) - Mkt
	max0<M_QP,data_t>(temp3, zk);				// zk
	vadd<M_QP,data_t>(temp2, Mtk, uk);			// uk
	}
	return;
}

void zk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t temp[M_QP], temp1[M_QP], temp2[M_QP];
    {
#pragma HLS LOOP_MERGE
    //  zk = max{0, c - uk - M*tk};
	vsub<M_QP,data_t>(c_qp, uk, temp);
	mvmult<M_QP,N_QP,data_t>(M_qp, tk, temp1);
	vsub<M_QP,data_t>(temp, temp1, temp2);
	max0<M_QP,data_t>(temp2, zk);
    }
	return;
}

void uk_update(data_t (&c_qp)[M_QP], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t temp[M_QP], temp1[M_QP], temp2[M_QP];
    {
#pragma HLS LOOP_MERGE
    // uk = uk + (M*tk + zk - c);
	vsub<M_QP,data_t>(zk, c_qp, temp);
	vadd<M_QP,data_t>(uk, temp, temp1);
	mvmult<M_QP,N_QP,data_t>(M_qp, tk, temp2);
	vadd<M_QP,data_t>(temp1, temp2, uk);
    }
	return;
}

