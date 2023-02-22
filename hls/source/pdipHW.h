#ifndef PDIPHW_H
#define PDIPHW_H
#include <iostream>

#include "specs.h"
#include "linearSolversHW.h"
#include "mmultHW.h"

using namespace std;

void pdipHW(elem H[N_HOR][N_HOR], elem h[N_HOR], elem M[I_AUX][N_HOR], elem c[I_AUX], elem tk[N_HOR], int iterPDIP, int iterLs, elem tol);
#endif
