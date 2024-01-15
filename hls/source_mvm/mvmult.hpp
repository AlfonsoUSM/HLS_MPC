
#ifndef MVMULT_H
#define MVMULT_H

#define N PARAM_N
#define M PARAM_M
typedef float T;

void mvmult_rowwise(T (&A)[N][M], T (&B)[M], T (&R)[N]);

void mvmult_columnwise(T (&A)[N][M], T (&B)[M], T (&R)[N]);


#endif // MVMULT_H
