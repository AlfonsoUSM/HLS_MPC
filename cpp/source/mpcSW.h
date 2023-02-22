#ifndef MPCSW_H
#define MPCSW_H

#include<iostream>

#include"specs.h"
#include"pdipSW.h"

using namespace std;

void mpcSW(elem x[N_SYS], elem yref[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[P_SYS], elem umax[P_SYS], elem Acal[N_SYS*N_HOR][N_SYS], elem AcalQOcal[N_SYS][N_HOR], elem H[N_HOR][N_HOR], elem M_hat[I_AUX][N_HOR], elem L_invLast[N_SYS+P_SYS], int iterPDIP, int iterLS, elem tol, elem u_star[M_SYS]);

void get_h_tilde(elem AcalQOcal[N_SYS][N_HOR], elem x_tilde[N_SYS], elem h_tilde[N_HOR]);
void get_c_hat(elem Acal[N_SYS*N_HOR][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_star[N_SYS], elem x_tilde[N_SYS], elem a_tilde[N_HOR], elem b_tilde[N_HOR], elem c_hat[I_AUX]);

#endif
