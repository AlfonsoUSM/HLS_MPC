
#ifndef MPC_H
#define MPC_H

#include "system.hpp"

/*!
@brief
*/
void mpc_sparse_iteration(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&x1)[N_SYS], data_t(&u0)[M_SYS]);

/*!
@brief
*/
void mpc_sparse_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&c_hat)[M_QP]);


#endif // MPC_H

