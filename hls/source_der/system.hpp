#ifndef SYSTEM_H
#define SYSTEM_H


////////// System //////////
#define M_SYS 2		// control inputs
#define N_SYS 4		// system states
#define P_SYS 2		// system references
#define A_SYS 10	// input constraints pairs
#define B_SYS 10	// state constraints pairs
#define C_SYS 2*(A_SYS+B_SYS)
#define D_SYS 2		// disturbances

#define N_HOR	2	// prediction horizon
#define DENSE//SPARSE//
typedef float data_t;
extern data_t Bpd[N_SYS][D_SYS];
extern data_t umax[A_SYS];
extern data_t umin[A_SYS];
extern data_t xmax[B_SYS];
extern data_t xmin[B_SYS];
extern data_t H[A_SYS][M_SYS];
extern data_t I[B_SYS][N_SYS];

////////// MPC Formulation //////////

#if defined DENSE		// Dense formulation
#define N_QP M_SYS*N_HOR	// Number of optimization decision variables
#define M_QP C_SYS*N_HOR	// Number of optimization constraints
//extern data_t c[5*N_HOR];
//extern data_t d[5*N_HOR];
//extern data_t e[5*N_HOR];
//extern data_t f[5*N_HOR];
//extern data_t D[5*N_HOR][N_SYS];
extern data_t F[N_SYS][N_QP];
extern data_t T_inv[N_SYS+M_SYS][N_SYS+M_SYS];
//extern data_t V[A_SYS*N_HOR][N_QP];
//extern data_t W[B_SYS*N_HOR][N_QP];
extern data_t WD[B_SYS*N_HOR][N_SYS];

#elif defined SPARSE	// Sparse formulation
#define N_QP ((N_HOR + 1)*N_SYS + N_HOR*M_SYS)					// Number of optimization variables
#define M_QP (2*(N_HOR + 1)*N_SYS + 2*N_HOR*A_SYS + 2*(N_HOR + 1)*B_SYS)	// Number of optimization constraints
extern data_t g[(2*N_QP)];
// f is (N_HOR+1)*N_SYSx1

#else
#error Something is wrong with how elem is defined
#endif

// Q[N_QP][N_QP] is not required
// q[N_QP] is not constant
extern data_t G[M_QP][N_QP];		// G_hat
// g[M_QP1] is not constant


////////// ADMM //////////

//rho is not required
extern data_t R_inv[N_QP][N_QP];
extern data_t P[N_QP][M_QP];		// RhoGt_neg


#endif // SYSTEM_H



