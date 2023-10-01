
#include "system.hpp"

// HOR = 5
#if defined DENSE


data_t H[M_QP][N_QP] =	{{3.853948414325714e-02, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {1.935940963448957e-05, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {3.747834637761116e-02, 3.853948414325714e-02, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {5.736585808335803e-05, 1.935940963448957e-05, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {3.644642233848572e-02, 3.747834637761116e-02, 3.853948414325714e-02, 0.000000000000000e+00, 0.000000000000000e+00},
					 {9.432584192836657e-05, 5.736585808335803e-05, 1.935940963448957e-05, 0.000000000000000e+00, 0.000000000000000e+00},
					 {3.544291108846664e-02, 3.644642233848572e-02, 3.747834637761116e-02, 3.853948414325714e-02, 0.000000000000000e+00},
					 {1.302681775996462e-04, 9.432584192836657e-05, 5.736585808335803e-05, 1.935940963448957e-05, 0.000000000000000e+00},
					 {3.446703404188156e-02, 3.544291108846664e-02, 3.644642233848572e-02, 3.747834637761116e-02, 3.853948414325714e-02},
					 {1.652208738960326e-04, 1.302681775996462e-04, 9.432584192836657e-05, 5.736585808335803e-05, 1.935940963448957e-05},
					 {-3.853948414325714e-02, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-1.935940963448957e-05, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-3.747834637761116e-02, -3.853948414325714e-02, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-5.736585808335803e-05, -1.935940963448957e-05, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-3.644642233848572e-02, -3.747834637761116e-02, -3.853948414325714e-02, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-9.432584192836657e-05, -5.736585808335803e-05, -1.935940963448957e-05, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-3.544291108846664e-02, -3.644642233848572e-02, -3.747834637761116e-02, -3.853948414325714e-02, -0.000000000000000e+00},
					 {-1.302681775996462e-04, -9.432584192836657e-05, -5.736585808335803e-05, -1.935940963448957e-05, -0.000000000000000e+00},
					 {-3.446703404188156e-02, -3.544291108846664e-02, -3.644642233848572e-02, -3.747834637761116e-02, -3.853948414325714e-02},
					 {-1.652208738960326e-04, -1.302681775996462e-04, -9.432584192836657e-05, -5.736585808335803e-05, -1.935940963448957e-05},
					 {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00},
					 {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00},
					 {-1.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-0.000000000000000e+00, -1.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -1.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -1.000000000000000e+00, -0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -1.000000000000000e+00}};

data_t a_neg[N_QP] =	{3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00};

data_t b[N_QP] =	{3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00, 3.000000000000000e+00};

data_t d[N_SYS*N_HOR] =	{-5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00, -5.000000000000000e+00, -2.000000000000000e+00};

data_t e[N_SYS*N_HOR] =	{5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 2.000000000000000e+00};


data_t D[N_SYS*N_HOR][N_SYS] =	{{9.724661707878113e-01, 0.000000000000000e+00},
					 {9.861689759418368e-04, 1.000000000000000e+00},
					 {9.456904530525208e-01, 0.000000000000000e+00},
					 {1.945185009390116e-03, 1.000000000000000e+00},
					 {9.196519255638123e-01, 0.000000000000000e+00},
					 {2.877795603126287e-03, 1.000000000000000e+00},
					 {8.943303823471069e-01, 0.000000000000000e+00},
					 {3.784727770835161e-03, 1.000000000000000e+00},
					 {8.697060346603394e-01, 0.000000000000000e+00},
					 {4.666688852012157e-03, 1.000000000000000e+00}};


data_t G[N_SYS][N_QP] =	{{2.151192724704742e-02, 2.146639861166477e-02, 2.141969650983810e-02, 2.137186191976070e-02, 2.132293581962585e-02},
					 {6.442398428916931e-01, 6.416124701499939e-01, 6.389900445938110e-01, 6.363726854324341e-01, 6.337603330612183e-01}};


data_t R_inv[N_QP][N_QP] =	{{2.477079391479492e+00, -1.164435222744942e-02, -9.981122799217701e-03, -8.329093456268311e-03, -6.685729604214430e-03},
					 {-1.164435315877199e-02, 2.478572130203247e+00, -1.011718809604645e-02, -8.417451754212379e-03, -6.726513616740704e-03},
					 {-9.981121867895126e-03, -1.011718716472387e-02, 2.480135917663574e+00, -8.514588698744774e-03, -6.773468106985092e-03},
					 {-8.329093456268311e-03, -8.417452685534954e-03, -8.514588698744774e-03, 2.481780052185059e+00, -6.826665252447128e-03},
					 {-6.685728672891855e-03, -6.726513616740704e-03, -6.773468106985092e-03, -6.826665718108416e-03, 2.483514308929443e+00}};


data_t W[N_QP][M_QP] =	{{-3.883658908307552e-03, -1.950865453181905e-06, -3.776727244257927e-03, -5.780809715361102e-06, -3.672739258036017e-03, -9.505301022727508e-06, -3.571614623069763e-03, -1.312724361923756e-05, -3.473274409770966e-03, -1.664945921220351e-05, 3.883658908307552e-03, 1.950865453181905e-06, 3.776727244257927e-03, 5.780809715361102e-06, 3.672739258036017e-03, 9.505301022727508e-06, 3.571614623069763e-03, 1.312724361923756e-05, 3.473274409770966e-03, 1.664945921220351e-05, -1.007709130644798e-01, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, 1.007709130644798e-01, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -3.883658908307552e-03, -1.950865453181905e-06, -3.776727244257927e-03, -5.780809715361102e-06, -3.672739258036017e-03, -9.505301022727508e-06, -3.571614623069763e-03, -1.312724361923756e-05, 0.000000000000000e+00, 0.000000000000000e+00, 3.883658908307552e-03, 1.950865453181905e-06, 3.776727244257927e-03, 5.780809715361102e-06, 3.672739258036017e-03, 9.505301022727508e-06, 3.571614623069763e-03, 1.312724361923756e-05, -0.000000000000000e+00, -1.007709130644798e-01, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, 0.000000000000000e+00, 1.007709130644798e-01, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -3.883658908307552e-03, -1.950865453181905e-06, -3.776727244257927e-03, -5.780809715361102e-06, -3.672739258036017e-03, -9.505301022727508e-06, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 3.883658908307552e-03, 1.950865453181905e-06, 3.776727244257927e-03, 5.780809715361102e-06, 3.672739258036017e-03, 9.505301022727508e-06, -0.000000000000000e+00, -0.000000000000000e+00, -1.007709130644798e-01, -0.000000000000000e+00, -0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.007709130644798e-01, 0.000000000000000e+00, 0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -3.883658908307552e-03, -1.950865453181905e-06, -3.776727244257927e-03, -5.780809715361102e-06, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 3.883658908307552e-03, 1.950865453181905e-06, 3.776727244257927e-03, 5.780809715361102e-06, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -1.007709130644798e-01, -0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.007709130644798e-01, 0.000000000000000e+00},
					 {-0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -3.883658908307552e-03, -1.950865453181905e-06, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 3.883658908307552e-03, 1.950865453181905e-06, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, -1.007709130644798e-01, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.007709130644798e-01}};


#else

// SPARSE

#endif