#ifndef MPCHW_H
#define MPCHW_H

#include<iostream>

#include"specs.h"
#include"pdipHW.h"

using namespace std;

void mpcHW(elem x_axi[N_SYS], elem r_axi[P_SYS], elem xmin_axi[N_SYS], elem xmax_axi[N_SYS], elem umin_axi[M_SYS], elem umax_axi[M_SYS], elem Acal_axi[N_SYS*N_HOR][N_SYS], elem AcalQOcal_axi[N_SYS][N_HOR], elem H_axi[N_HOR][N_HOR], elem M_hat_axi[I_AUX][N_HOR], elem L_invLast_axi[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u_axi[M_SYS]);

void get_h_tilde(elem AcalQOcal[N_SYS][N_HOR], elem x_tilde[N_SYS], elem h_tilde[N_HOR]);
void get_c_hat(elem Acal[N_SYS*N_HOR][N_SYS], elem xmax[N_SYS], elem xmin[N_SYS], elem x_star[N_SYS], elem x_tilde[N_SYS], elem a[N_HOR], elem b[N_HOR], elem c_hat[6*N_HOR]);

#endif
