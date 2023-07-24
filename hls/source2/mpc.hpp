
#ifndef MPC_H
#define MPC_H

#include "system.hpp"

/*!
@brief
*/
void mpc_sparse_admm_iteration(const data_t (&r0)[P_SYS], const data_t (&x0)[N_SYS], const data_t (&M_hat)[M_QP][N_QP], const data_t (&RhoMt_neg)[N_QP][M_QP], const data_t (&h_qp)[N_QP], const data_t (&R_inv)[N_QP][N_QP], data_t (&x1)[N_SYS], data_t (&u0)[M_SYS]);

/*!
@brief
*/
void mpc_sparse_constraint(const data_t (&r0)[P_SYS], const data_t (&x0)[N_SYS], data_t (&c_hat)[M_QP]);


#endif // MPC_H

