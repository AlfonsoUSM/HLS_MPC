#ifndef SYSTEM_H
#define SYSTEM_H

typedef float data_t;
#define SPARSE

#define N_HOR	4//16		// prediction horizon
#define N_IT	10		// number of iterations of QP solver

////////// System //////////
#define N_SYS 2		// system states
#define M_SYS 1		// control inputs
#define P_SYS 1		// system outputs

//extern data_t d[2*N_SYS];		// [xmax; -xmin]
//extern data_t e[2*M_SYS];		// [umax; -umin]

/*
 * A is NxN
 * B is NxM
 * C is PxN
 *
 * Gamma is MxM
 * Omega is NxN
 */

////////// MPC Formulation //////////
//#if !defined DENSE and !defined SPARSE
//#define DENSE
//#endif

#if defined DENSE		// Dense formulation
#define N_QP N_SYS*N_HOR							// Number of optimization values
#define M_QP (2 * N_HOR * (N_SYS + M_SYS))			// Number of optimization constraints
data_t a_tilde[N_HOR];
data_t b_tilde[N_HOR];

#elif defined SPARSE	// Sparse formulation
#define N_QP (N_HOR * (N_SYS + M_SYS) + N_SYS)					// Number of optimization variables
#define M_QP (4 * ( N_HOR + 1) * N_SYS + 2 * N_HOR * M_SYS)		// Number of optimization constraints
extern data_t g[(2*N_QP)];
// f is (N_HOR+1)*N_SYSx1

#else
#error Something is wrong with how elem is defined
#endif

// Q[N_QP][N_QP] is not required
extern data_t q[N_QP];				// q
extern data_t H[M_QP][N_QP];		// H_hat
// h[M_QP1] is not constant


////////// ADMM //////////

//rho is not required
extern data_t R_inv[N_QP][N_QP];
extern data_t RhoHt_neg[N_QP][M_QP];


#endif // SYSTEM_H



