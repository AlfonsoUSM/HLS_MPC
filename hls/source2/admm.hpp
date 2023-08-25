#ifndef ADMM_H
#define ADMM_H

#include "system.hpp"
#include "utils.hpp"

extern data_t Rho;

extern data_t tk_admm[N_QP];
extern data_t zk_admm[M_QP];
extern data_t uk_admm[M_QP];

/*!
@brief 	Alternating Method of Multipliers for solving quadratic convex optimization.
        Quadratic programing (QP) problem solution:
                    min 1/2 tk'*H*tk + h'*tk
                    s.t     M*tk <= c
        Arrange the constraints into equality constraints for ADMM formulation, being finally:
                    min 1/2 tk'*H*tk + h'*tk + g(zk)
                    s.t     M*tk + zk = c
                            zk >= 0
        where g(zk) is the indicator function of Z
            g(zk) = 0    if zk in Z
            g(zk) = ∞k    if any component of zk not in Z

@tparam N_QP	Number of optimization values
@tparam M_QP	Number of optimization constraints
@param  c_qp	M_QPx1 vector with constraints
@return uk		N_QPx1 optimal solutions vector
*/
void qp_admm(const data_t (&c_qp)[M_QP], const data_t (&C_qp)[M_QP][N_QP], const data_t (&RhoMt_neg)[N_QP][M_QP],const  data_t (&h_qp)[N_QP], const data_t (&R_inv)[N_QP][N_QP]);

/*!
@brief
*/
void tk_update(const data_t (&c_qp)[M_QP], const data_t (&RhoMt_neg)[N_QP][M_QP], const data_t (&h_qp)[N_QP], const data_t (&R_inv)[N_QP][N_QP], data_t (&tk)[N_QP], const data_t (&zk)[M_QP], const data_t (&uk)[M_QP]);

/*!
@brief
*/
void zk_uk_update (const data_t (&c_qp)[M_QP], const data_t (&C_qp)[M_QP][N_QP], const data_t (&tk)[N_QP], data_t (&zk)[M_QP], data_t (&uk)[M_QP]);

/*!
@brief
*/
void zk_update(const data_t (&c_qp)[M_QP], const data_t (&C_qp)[M_QP][N_QP], const data_t (&tk)[N_QP], data_t (&zk)[M_QP], const data_t (&uk)[M_QP]);

/*!
@brief
*/
void uk_update(const data_t (&c_hat)[M_QP], const data_t (&C_qp)[M_QP][N_QP], const data_t (&tk)[N_QP], const data_t (&zk)[M_QP], data_t (&uk)[M_QP]);



#endif // ADMM_H

