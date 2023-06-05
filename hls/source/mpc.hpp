
#ifndef MPC_H
#define MPC_H

#include "system.hpp"

/*!
@brief
*/
void mpc_sparse_admm_iteration(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS],data_t (&theta)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]);

/*!
@brief
*/
void mpc_sparse_constraint(data_t (&r0)[P_SYS], data_t (&x0)[N_SYS], data_t (&c_hat)[M_QP]);


#endif // MPC_H

