
#include "mvmult.hpp"

void mvmult_rowwise(T (&A)[N][M], T (&B)[M], T (&R)[N]){
	outer_loop: for(int i = 0; i < N; ++i){
		inner_loop: for(int j = 0; j < M; ++j){
			if (j == 0){
				R[i] = A[i][j] * B[j];
			}
			else{
				R[i] += A[i][j] * B[j];
			}
		}
	}
	return;
}

void mvmult_columnwise(T (&A)[N][M], T (&B)[M], T (&R)[N]){
	outer_loop: for(int j = 0; j < M; ++j){
		inner_loop: for(int i = 0; i < N; ++i){
			if (j == 0){
				R[i] = A[i][j] * B[j];
			}
			else{
				R[i] += A[i][j] * B[j];
			}
		}
	}
	return;
}
