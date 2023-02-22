#ifndef MMULTHW_H
#define MMULTHW_H

#include "specs.h"
#define MAX_SIZE I_AUX

elem adderTree (elem input[MAX_SIZE]);
void mmultHW (elem A[MAX_SIZE][MAX_SIZE], elem B[MAX_SIZE][MAX_SIZE], elem C[MAX_SIZE][MAX_SIZE], int m, int n, int p);
void mmult_N_I_N (elem A[N_HOR][I_AUX], elem B[I_AUX][N_HOR], elem C[N_HOR][N_HOR]);
void mmult_N_N_1 (elem A[N_HOR][N_HOR], elem B[N_HOR], elem C[N_HOR]);
void mmult_N_I_1 (elem A[N_HOR][I_AUX], elem B[I_AUX], elem C[N_HOR]);
void mmult_I_N_1 (elem A[I_AUX][N_HOR], elem B[N_HOR], elem C[I_AUX]);

#endif
