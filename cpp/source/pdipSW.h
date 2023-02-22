#ifndef PDIPSW_H
#define PDIPSW_H
#include <iostream>

#include "specs.h"
#include "linearSolversSW.h"
#include "mmultSW.h"

using namespace std;


void pdipSW(elem H[N_HOR][N_HOR], elem h[N_HOR], elem M[I_AUX][N_HOR], elem c[I_AUX], elem tk[N_HOR], int iterPDIP, int iterLs, elem tol);
#endif
