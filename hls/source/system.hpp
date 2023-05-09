#ifndef SYSTEM_H
#define SYSTEM_H

//#include "specs.h"
//typedef int data_t;
typedef float data_t;

#define N_HOR	4		// prediction horizon
#define N_IT	10		// number of iterations of QP solver

////////// System //////////
#define N_SYS 2		// system states
#define M_SYS 1		// control inputs
#define P_SYS 1		// system outputs

/*
 * A is NxN
 * B is NxM
 * C is PxN
 *
 * Gamma is MxM
 * Omega is NxN
 */

////////// Formulation //////////
#if !defined DENSE and !defined SPARSE
#define SPARSE
#endif

#if defined DENSE // Dense formulation
#define N_QP N_SYS*N_HOR							// Number of optimization values
data_t a_tilde[N_HOR];
data_t b_tilde[N_HOR];

#elif defined SPARSE // Sparse formulation
#define N_QP (N_HOR * (N_SYS + M_SYS) + N_SYS)					// Number of optimization variables
#define M_QP (4 * ( N_HOR + 1) * N_SYS + 2 * N_HOR * M_SYS)	// Number of optimization constraints
extern data_t g[(2*N_QP)];
// f is (N_HOR+1)*N_SYSx1

#else
#error Something is wrong with how elem is defined
#endif

extern data_t H_qp[N_QP][N_QP];
extern data_t h_qp[N_QP];
extern data_t M_qp[M_QP][N_QP];
// c_qp is M_QPx1

extern data_t Rho;
extern data_t R_inv[N_QP][N_QP];
extern data_t RhoMt_neg[N_QP][M_QP];

#endif // SYSTEM_H


