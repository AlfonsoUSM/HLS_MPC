#ifndef LINEARSOLVERSHW_H
#define LINEARSOLVERSHW_H

#include "specs.h"
#include <stdio.h>
#include <math.h>

// linear solvers
void minresHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zko[N_HOR], elem zk[N_HOR], elem iter, elem tol);

void cgradHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zko[N_HOR], elem zk[N_HOR], elem tol);

void cholHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zk[N_HOR]);

#endif
