
#ifndef MPC_H
#define MPC_H

#include "system.hpp"

////// MPC //////

/*!
@brief
*/
void mpc(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&d0)[D_SYS], data_t (&u0)[M_SYS], int IT);

/*!
@brief
*/
void mpc_dense_constraint(data_t (&x0)[N_SYS], data_t (&r0)[P_SYS], data_t (&d0)[D_SYS], data_t (&inf)[N_SYS+M_SYS], data_t (&q)[N_QP], data_t (&g)[M_QP]);

/*!
@brief
*/
void mpc_sparse_constraint(data_t (&x0)[N_SYS], data_t (&g)[M_QP]);

////// ADMM //////

extern data_t tk_admm[N_QP];
extern data_t zk_admm[M_QP];
extern data_t uk_admm[M_QP];

/*!
@brief 	Alternating Method of Multipliers for solving quadratic convex optimization.
        Quadratic programing (QP) problem solution:
                    min 1/2 tk'*Q*tk + q'*tk
                    s.t     H*tk <= h
        Arrange the constraints into equality constraints for ADMM formulation, being finally:
                    min 1/2 tk'*Q*tk + q'*tk + g(zk)
                    s.t     H*tk + zk = h
                            zk >= 0
        where g(zk) is the indicator function of Z
            g(zk) = 0    if zk in Z
            g(zk) = ∞k    if any component of zk not in Z

@tparam N_QP	Number of optimization values
@tparam M_QP	Number of optimization constraints
@param  h		M_QPx1 vector with constraints
@return uk		N_QPx1 optimal solutions vector
*/
void qp_admm(data_t (&q)[N_QP], data_t (&h)[M_QP], int IT);

#endif // MPC_H

