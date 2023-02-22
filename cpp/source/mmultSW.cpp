#include "mmultSW.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>


void mmultSW (elem *A, elem *B, elem *C, int m, int n, int p){

    for (int i = 0 ; i < m ; i++){
        for(int j = 0 ; j < p ; j++){
        	elem sum = 0;
            for(int k = 0; k < n; k++){
            	sum += A[i*n +k]*B[k*p +j];
            }
            C[i*p +j]= sum;
        }
    }

}
