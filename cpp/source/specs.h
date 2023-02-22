#ifndef SPECS_H
#define SPECS_H

#ifndef N_HOR
#define N_HOR 4
#endif

#define N_SYS 2
#define M_SYS 1
#define P_SYS 1
#define I_AUX (N_SYS*N_HOR*2+2*N_HOR)

// type of data (only one at a time)
#if !defined ELEM_DOULE and !defined ELEM_FLOAT
//#define ELEM_DOUBLE
#define ELEM_FLOAT
#endif

#if defined ELEM_DOUBLE
typedef double elem;
#elif defined ELEM_FLOAT
typedef float elem;
#else
#error Something is wrong with how elem is defined
#endif

// linear solver
#ifndef LS
#define LS 3
#endif

#define MINRES  1
#define CGRAD   2
#define CHOL  3

#endif