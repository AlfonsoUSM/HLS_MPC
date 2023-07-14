#include "system.hpp"


data_t d[2*N_SYS] = {5, 2, 5, 2};	// [xmax; -xmin]
data_t e[2*M_SYS] = {10, 10};		// [umax; -umin]

data_t H_qp[N_QP][N_QP] = {2};
data_t h_qp[N_QP] = {2};
data_t C_qp[M_QP][N_QP] = {2};

// sparse formulation
// g = [d; d; d... dN; e; e; e...]
data_t g[(2*N_QP)] = {2};

