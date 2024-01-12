#ifndef SYSTEM_H
#define SYSTEM_H


////////// System //////////
#define N_SYS 2		// system states
#define M_SYS 1		// control inputs
#define P_SYS 1		// system outputs

#define N_HOR	4//HOR_SIZE//16//		// prediction horizon
#define DENSE//SPARSE//
typedef float data_t;

////////// MPC Formulation //////////

#if defined DENSE		// Dense formulation
#define N_QP M_SYS*N_HOR							// Number of optimization values
#define M_QP (2 * N_HOR * (N_SYS + M_SYS))			// Number of optimization constraints
//extern data_t c[N_QP];
//extern data_t d_neg[N_QP];
//extern data_t e[N_SYS*N_HOR];
//extern data_t f[N_SYS*N_HOR];
extern data_t D[N_SYS*N_HOR][N_SYS];
extern data_t F[N_SYS][N_QP];
extern data_t T_inv[N_SYS+M_SYS][N_SYS+M_SYS];
extern data_t xmin[N_SYS];
extern data_t xmax[N_SYS];
extern data_t umin[M_SYS];
extern data_t umax[M_SYS];

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
extern data_t G[M_QP][N_QP];		// G_hat
// g[M_QP1] is not constant


////////// ADMM //////////

//rho is not required
extern data_t R_inv[N_QP][N_QP];
extern data_t P[N_QP][M_QP];		// RhoHt_neg


#endif // SYSTEM_H



