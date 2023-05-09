#include "system.hpp"

data_t H_qp[N_QP][N_QP] = {2};
data_t h_qp[N_QP] = {2};
data_t M_qp[M_QP][N_QP] = {2};
data_t g[(2*N_QP)] = {2};

// admm
data_t Rho = 2;
data_t R_inv[N_QP][N_QP] = {2};
data_t RhoMt_neg[N_QP][M_QP] = {2};
