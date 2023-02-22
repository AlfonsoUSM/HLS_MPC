#include "mmultHW.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>


elem adderTree (elem input[MAX_SIZE]){
#pragma HLS INLINE
//#pragma HLS unroll
	for(int i=1; i<MAX_SIZE; i = i*2){
		for(int j=0; j+i<(MAX_SIZE); j=(j+i*2)){
			input[j] = input[j] + input[j+i];
		}
	}
	return input[0];
}

void mmultHW (elem A[MAX_SIZE][MAX_SIZE], elem B[MAX_SIZE][MAX_SIZE], elem C[MAX_SIZE][MAX_SIZE], int m, int n, int p){
//#pragma HLS ARRAY_PARTITION variable=A cyclic factor=24 dim=1

    //loop tripcount constant
    const int c_size = MAX_SIZE;

    /*
    const int lim_loop_r = (m<n)? n: m;
    const int lim_loop_c = (n<p)? p: n;

    //Local memory to store input matrices
    float local_in1[MAX_SIZE][MAX_SIZE];
    float local_in2[MAX_SIZE][MAX_SIZE];
    //float local_out[MAX_SIZE][MAX_SIZE];


    //Physical implementation of memories have only a limited number of read/write
    //ports, that can be overcome by using the ARRAY_PARTITION pragma
    #pragma HLS ARRAY_PARTITION variable=local_in1 complete dim=2
    #pragma HLS ARRAY_PARTITION variable=local_in2 complete dim=1


//LOOP MAX INPUT MATRIX

    for (int r=0; r<lim_loop_r; r++){
	#pragma HLS LOOP_TRIPCOUNT min=c_size max=c_size
    	for (int c=0; c<lim_loop_c; c++){
	#pragma HLS LOOP_TRIPCOUNT min=c_size max=c_size
	#pragma HLS PIPELINE

    		if ((r<m) & (c<n)) {
    			local_in1[r][c] = A[r*n+c];
    		}
    		//else{
    		//	local_in1[r][c] = 0;
    		//}
    		//load b
    		if ((r<n) & (c<p)) {
    			local_in2[r][c] = B[r*p+c];
    		}
    		//else{
    		//	local_in2[r][c] = 0;
    		//}
    	}
    }

*/

    //Reads the input_data from local memory, performs the
    //computations and writes the data to local memory.
    mult_1: for (int i = 0 ; i < m ; i++){
	#pragma HLS LOOP_TRIPCOUNT min=c_size max=c_size
        mult_2: for(int j = 0 ; j < p ; j++){
		#pragma HLS LOOP_TRIPCOUNT min=c_size max=c_size
        //Pipelining a loop results in automatic unrolling of inner loops by the HLS compiler.
		#pragma HLS PIPELINE II=1
        	elem multArray[MAX_SIZE];

            mult_3: for(int k = 0; k < MAX_SIZE; k++){
            	if (k<n){
            		//multArray[k] = local_in1[i][k]*local_in2[k][j];
            		multArray[k] = A[i][k]*B[k][j];
            	}
            	else{
            		//multArray[k] = 0;
            		multArray[k] = 0;
            	}

            }

            //C[i*p+j] = adderTree (multArray);
            C[i][j] = adderTree (multArray);
        }
    }
    /*
    //Burst write from output matrix local_out to DDR memory
    write_out: for(int iter = 0, i = 0, j = 0; iter < ROWS * COLS; iter++, j++){
    //Partitioning the output matrix local_out is not useful because the same element is
    //accessed and processed in all iterations of the innermost loop.
    #pragma HLS PIPELINE
    #pragma HLS LOOP_TRIPCOUNT min=c_size*c_size max=c_size*c_size
        if(j == ROWS){ j = 0; i++; }
        C[iter] = local_out[i][j];
    }
	*/

}


void mmult_N_I_N (elem A[N_HOR][I_AUX], elem B[I_AUX][N_HOR], elem C[N_HOR][N_HOR]){
#pragma HLS INLINE
//#pragma HLS ARRAY_PARTITION variable=A complete dim=2
//#pragma HLS ARRAY_PARTITION variable=B complete dim=1
    mult_1: for (int r = 0 ; r < N_HOR ; r++){ 		// rows of A
		for(int c = 0 ; c < N_HOR ; c++){			// cols of B
		#pragma HLS PIPELINE II=1
			elem multArray[MAX_SIZE];
			for(int k = 0; k < MAX_SIZE; k++){
				if (k<I_AUX){
					multArray[k] = A[r][k]*B[k][c];
				}
				else{
					multArray[k] = 0;
				}
			}
			C[r][c] = adderTree (multArray);
		}
	}
}

void mmult_N_N_1 (elem A[N_HOR][N_HOR], elem B[N_HOR], elem C[N_HOR]){
#pragma HLS INLINE
//#pragma HLS ARRAY_PARTITION variable=A complete dim=2
//#pragma HLS ARRAY_PARTITION variable=B complete dim=1
    mult_1: for (int r = 0 ; r < N_HOR ; r++){ 		// rows of A
		#pragma HLS PIPELINE II=1
    	elem multArray[MAX_SIZE];
		for(int k = 0; k < MAX_SIZE; k++){
			if (k<N_HOR){
				multArray[k] = A[r][k]*B[k];
			}
			else{
				multArray[k] = 0;
			}
		}
		C[r] = adderTree (multArray);
	}
}


void mmult_N_I_1 (elem A[N_HOR][I_AUX], elem B[I_AUX], elem C[N_HOR]){
#pragma HLS INLINE
//#pragma HLS ARRAY_PARTITION variable=A complete dim=2
//#pragma HLS ARRAY_PARTITION variable=B complete dim=1
    mult_1: for (int r = 0 ; r < N_HOR ; r++){ 		// rows of A
		#pragma HLS PIPELINE II=1
    	elem multArray[MAX_SIZE];
		for(int k = 0; k < MAX_SIZE; k++){
			if (k<I_AUX){
				multArray[k] = A[r][k]*B[k];
			}
			else{
				multArray[k] = 0;
			}
		}
		C[r] = adderTree (multArray);
	}
}

void mmult_I_N_1 (elem A[I_AUX][N_HOR], elem B[N_HOR], elem C[I_AUX]){
#pragma HLS INLINE
//#pragma HLS ARRAY_PARTITION variable=A complete dim=2
//#pragma HLS ARRAY_PARTITION variable=B complete dim=1
    mult_1: for (int r = 0 ; r < I_AUX ; r++){ 		// rows of A
		#pragma HLS PIPELINE II=1
    	elem multArray[MAX_SIZE];
		for(int k = 0; k < MAX_SIZE; k++){
			if (k<N_HOR){
				multArray[k] = A[r][k]*B[k];
			}
			else{
				multArray[k] = 0;
			}
		}
		C[r] = adderTree (multArray);
	}
}
