#ifndef SYSTEM_H
#define SYSTEM_H


////////// System //////////
#define N_SYS 2		// system states
#define M_SYS 1		// control inputs
#define P_SYS 1		// system outputs

#define N_HOR	HOR_SIZE//4//16		// prediction horizon
#define DENSE//SPARSE//
typedef float data_t;

////////// MPC Formulation //////////

#if defined DENSE		// Dense formulation
#define N_QP M_SYS*N_HOR							// Number of optimization values
#define M_QP (2 * N_HOR * (N_SYS + M_SYS))			// Number of optimization constraints
extern data_t a_neg[N_QP];
extern data_t b[N_QP];
extern data_t d[N_SYS*N_HOR];
extern data_t e[N_SYS*N_HOR];
extern data_t D[N_SYS*N_HOR][N_SYS];
extern data_t G[N_SYS][N_QP];

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
extern data_t W[N_QP][M_QP];		// RhoHt_neg


#endif // SYSTEM_H



