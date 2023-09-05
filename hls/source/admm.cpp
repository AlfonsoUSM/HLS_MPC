
#include "admm.hpp"
#include "utils.hpp"

data_t Rho = 62.9634;

data_t R_inv[N_QP][N_QP] = {{ 3.0803e-03,-8.6322e-07, 1.1617e-03, 4.3977e-07, 4.4042e-04, 2.9486e-07, 1.7305e-04,-2.6384e-07,-8.4035e-05, 1.3520e-06,-7.0511e-05,-2.6505e-05,-9.8146e-06,-3.2397e-06},
         	 	 	 	 	{-8.6322e-07, 3.0221e-03,-7.1661e-07, 1.1491e-03,-5.8172e-07, 4.3430e-04,-7.6465e-07, 1.5726e-04, 1.6151e-06,-3.8731e-05,-3.1454e-08,-9.3781e-09,-1.3001e-08,-3.5802e-08},
 							{ 1.1617e-03,-7.1661e-07, 3.5225e-03,-1.0081e-06, 1.3354e-03, 9.0587e-09, 5.2472e-04,-1.1204e-06,-2.5481e-04, 4.1784e-06, 9.2005e-05,-8.0366e-05,-2.9759e-05,-9.8232e-06},
 							{ 4.3977e-07, 1.1491e-03,-1.0081e-06, 3.4564e-03,-1.3148e-06, 1.3063e-03,-2.1291e-06, 4.7303e-04, 4.7753e-06,-1.1650e-04,-1.0642e-08,-5.4387e-08,-4.8800e-08,-1.1089e-07},
 							{ 4.4042e-04,-5.8172e-07, 1.3354e-03,-1.3148e-06, 3.6067e-03,-2.4893e-06, 1.4172e-03,-3.9364e-06,-6.8820e-04, 1.1509e-05, 3.4881e-05, 8.8750e-05,-8.0374e-05,-2.6530e-05},
 							{ 2.9486e-07, 4.3430e-04, 9.0587e-09, 1.3063e-03,-2.4893e-06, 3.4951e-03,-5.2923e-06, 1.2656e-03, 1.2580e-05,-3.1169e-04, 6.1665e-09,-5.3781e-08,-1.5348e-07,-3.0425e-07},
 							{ 1.7305e-04,-7.6465e-07, 5.2472e-04,-2.1291e-06, 1.4172e-03,-5.2923e-06, 3.7702e-03,-1.2909e-05,-1.8309e-03, 3.1219e-05, 1.3706e-05, 3.4872e-05, 9.1979e-05,-7.0580e-05},
 							{-2.6384e-07, 1.5726e-04,-1.1204e-06, 4.7303e-04,-3.9364e-06, 1.2656e-03,-1.2909e-05, 3.3338e-03, 3.2637e-05,-8.2104e-04,-2.7118e-08,-9.4155e-08,-3.0922e-07,-8.2077e-07},
 							{-8.4035e-05, 1.6151e-06,-2.5481e-04, 4.7753e-06,-6.8820e-04, 1.2580e-05,-1.8309e-03, 3.2637e-05, 4.8606e-03,-8.4461e-05,-6.6555e-06,-1.6934e-05,-4.4667e-05,-1.1843e-04},
 							{ 1.3520e-06,-3.8731e-05, 4.1784e-06,-1.1650e-04, 1.1509e-05,-3.1169e-04, 3.1219e-05,-8.2104e-04,-8.4461e-05, 2.1579e-03, 1.0861e-07, 2.8254e-07, 7.6022e-07, 2.0545e-06},
 							{-7.0511e-05,-3.1454e-08, 9.2005e-05,-1.0642e-08, 3.4881e-05, 6.1665e-09, 1.3706e-05,-2.7118e-08,-6.6555e-06, 1.0861e-07, 7.9292e-03,-2.0991e-06,-7.7730e-07,-2.5658e-07},
 							{-2.6505e-05,-9.3781e-09,-8.0366e-05,-5.4387e-08, 8.8750e-05,-5.3781e-08, 3.4872e-05,-9.4155e-08,-1.6934e-05, 2.8254e-07,-2.0991e-06, 7.9295e-03,-1.9778e-06,-6.5283e-07},
 							{-9.8146e-06,-1.3001e-08,-2.9759e-05,-4.8800e-08,-8.0374e-05,-1.5348e-07, 9.1979e-05,-3.0922e-07,-4.4667e-05, 7.6022e-07,-7.7730e-07,-1.9778e-06, 7.9296e-03,-1.7219e-06},
 							{-3.2397e-06,-3.5802e-08,-9.8232e-06,-1.1089e-07,-2.6530e-05,-3.0425e-07,-7.0580e-05,-8.2077e-07,-1.1843e-04, 2.0545e-06,-2.5658e-07,-6.5283e-07,-1.7219e-06, 7.9303e-03}};

data_t RhoHt_neg[N_QP][M_QP] = {-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-61.23,-0.062,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 61.23, 0.062,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
								   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-61.23,-0.062,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 61.23, 0.062,   0.0,   0.0,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-61.23,-0.062,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 61.23, 0.062,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-61.23,-0.062,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 61.23, 0.062,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 62.96,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-2.427,-0.001,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 2.427, 0.001,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-2.427,-0.001,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 2.427, 0.001,   0.0,   0.0,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-2.427,-0.001,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 2.427, 0.001,   0.0,   0.0,
								   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-62.96, 62.96,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,-2.427,-0.001,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 2.427, 0.001};

data_t tk_admm[N_QP] = {0};
data_t zk_admm[M_QP] = {0};
data_t uk_admm[M_QP] = {0};

void qp_admm(data_t (&h)[M_QP], data_t (&H)[M_QP][N_QP], data_t (&q)[N_QP]){
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

/*
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
*/
