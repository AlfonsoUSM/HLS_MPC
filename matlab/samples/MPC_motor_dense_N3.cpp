
#include "system.hpp"

// HOR = 5
#if defined DENSE


data_t H[M_QP][N_QP] =	{{3.853948414325714e-02, 0.000000000000000e+00, 0.000000000000000e+00},
			 {1.935940963448957e-05, 0.000000000000000e+00, 0.000000000000000e+00},
			 {3.747834637761116e-02, 3.853948414325714e-02, 0.000000000000000e+00},
			 {5.736585808335803e-05, 1.935940963448957e-05, 0.000000000000000e+00},
			 {3.644642233848572e-02, 3.747834637761116e-02, 3.853948414325714e-02},
			 {9.432584192836657e-05, 5.736585808335803e-05, 1.935940963448957e-05},
			 {-3.853948414325714e-02, -0.000000000000000e+00, -0.000000000000000e+00},
			 {-1.935940963448957e-05, -0.000000000000000e+00, -0.000000000000000e+00},
			 {-3.747834637761116e-02, -3.853948414325714e-02, -0.000000000000000e+00},
			 {-5.736585808335803e-05, -1.935940963448957e-05, -0.000000000000000e+00},
			 {-3.644642233848572e-02, -3.747834637761116e-02, -3.853948414325714e-02},
			 {-9.432584192836657e-05, -5.736585808335803e-05, -1.935940963448957e-05},
			 {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
			 {0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00},
			 {0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00},
			 {-1.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
			 {-0.000000000000000e+00, -1.000000000000000e+00, -0.000000000000000e+00},
			 {-0.000000000000000e+00, -0.000000000000000e+00, -1.000000000000000e+00}};

data_t a_neg[N_QP] =	{3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00};

data_t b[N_QP] =	{3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00};

data_t d[N_SYS*N_HOR] =	{-5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00};

data_t e[N_SYS*N_HOR] =	{5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00};


data_t D[N_SYS*N_HOR][N_SYS] =	{{9.724661707878113e-01, 0.000000000000000e+00},
			 {9.861689759418368e-04, 1.000000000000000e+00},
			 {9.456904530525208e-01, 0.000000000000000e+00},
			 {1.945185009390116e-03, 1.000000000000000e+00},
			 {9.196519255638123e-01, 0.000000000000000e+00},
			 {2.877795603126287e-03, 1.000000000000000e+00}};


data_t G[N_SYS][N_QP] =	{{2.133546210825443e-02, 2.129026316106319e-02, 2.124390192329884e-02},
			 {6.389900445938110e-01, 6.363726854324341e-01, 6.337603330612183e-01}};


data_t R_inv[N_QP][N_QP] =	{{2.481985330581665e+00, -8.590785786509514e-03, -6.836050655692816e-03},
			 {-8.590785786509514e-03, 2.483654022216797e+00, -6.880118977278471e-03},
			 {-6.836051587015390e-03, -6.880118977278471e-03, 2.485409021377563e+00}};


data_t W[N_QP][M_QP] =	{{-3.877639304846525e-03, -1.947841610672185e-06, -3.770873183384538e-03, -5.771849828306586e-06, -3.667046548798680e-03, -9.490568118053488e-06, 3.877639304846525e-03, 1.947841610672185e-06, 3.770873183384538e-03, 5.771849828306586e-06, 3.667046548798680e-03, 9.490568118053488e-06, -1.006147190928459e-01, -0.000000000000000e+00, -0.000000000000000e+00, 1.006147190928459e-01, 0.000000000000000e+00, 0.000000000000000e+00},
			 {-0.000000000000000e+00, -0.000000000000000e+00, -3.877639304846525e-03, -1.947841610672185e-06, -3.770873183384538e-03, -5.771849828306586e-06, 0.000000000000000e+00, 0.000000000000000e+00, 3.877639304846525e-03, 1.947841610672185e-06, 3.770873183384538e-03, 5.771849828306586e-06, -0.000000000000000e+00, -1.006147190928459e-01, -0.000000000000000e+00, 0.000000000000000e+00, 1.006147190928459e-01, 0.000000000000000e+00},
			 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -3.877639304846525e-03, -1.947841610672185e-06, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 3.877639304846525e-03, 1.947841610672185e-06, -0.000000000000000e+00, -0.000000000000000e+00, -1.006147190928459e-01, 0.000000000000000e+00, 0.000000000000000e+00, 1.006147190928459e-01}};


#else

// SPARSE

#endif