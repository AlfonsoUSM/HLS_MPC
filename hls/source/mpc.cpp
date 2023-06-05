
#include "mpc.hpp"
#include "admm.hpp"
#include "utils.hpp"

/*
 * reference r0 = input()
 * state0 x0 = sense()
 * (u0, x1) = mpc_iteration(r0, x0)
 * actuate(u0)
 */


void mpc_sparse_admm_iteration(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]){
	data_t c_hat[M_QP] = {0};	// inicializar a 0 bien
	// constraint c_hat = constraint(x0, r0)
	mpc_sparse_constraint(r0, x0, c_hat);
	// optimized theta = qp_solver(H,h_nau, C_hat, c_hat)
	qp_admm(c_hat, tk, zk, uk);
	return;
}


// currently the reference is 0
void mpc_sparse_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&c_hat)[M_QP]){
	// follow reference currently not implemented
	// c_hat = [g; f; -f], where f = [-x0; 0; 0; 0...] is (N_HOR+1)*N_SYS
	for (int i=0; i<(2*N_QP); i++){
		c_hat[i] = g[i];
	}
	for (int i=0; i<N_SYS; i++){
		c_hat[(2*N_QP + i)] = -x0[i];
	}
	for (int i=0; i<N_SYS; i++){
		c_hat[(2*N_QP + N_HOR*N_SYS + N_SYS + i)] = x0[i];
	}
	return;
}

/*
void mpc_dense_iteration(){
	// reference r = input()
	// state0 x0 = sense()
	// constraint c_hat = constraint(x0, r)
	// actuation u = qp_solver(H,h_nau, M_hat, c_hat)
	// estimate x1
	// actuate()

	return;
}
*/
