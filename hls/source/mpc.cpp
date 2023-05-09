
#include "mpc.hpp"
#include "admm.hpp"
#include "utils.hpp"

/*
 * reference r0 = input()
 * state0 x0 = sense()
 * (u0, x1) = mpc_iteration(r0, x0)
 * actuate(u0)
 */


void mpc_sparse_iteration(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&x1)[N_SYS], data_t(&u0)[M_SYS]){
	data_t c_hat[M] = {2};
	data_t theta[N];
	// constraint c_hat = constraint(x0, r0)
//	mpc_sparse_constraint(r0, x0, c_hat);
	// actuation theta = qp_solver(H,h_nau, M_hat, c_hat)
	qp_admm(c_hat, theta);
	// estimate x1
	// x1 =
	// u0 =
	return;
}

void mpc_sparse_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&c_hat)[M]){

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
