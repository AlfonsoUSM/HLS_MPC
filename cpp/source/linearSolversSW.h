#ifndef LINEARSOLVERSSW_H
#define LINEARSOLVERSSW_H

#include "specs.h"
#include <stdio.h>
#include <math.h>

// linear solvers
void minresSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem iter, elem tol);

void cgradSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem tol);

void cholSW(elem* Ak, elem* bk, elem* zk);

#endif
