#pragma once

#include "ADMM.hpp"
#include "ADMM_wrapper.hpp"

/*!
@brief  wrapper for ADMM solver

@tparam N   Number of optimization values
@tparam M   Number of systems constraints
@tparam T   Data type
@param  P   NxN Cost matrix
@param  q   Nx1 Cost vector
@param  A  MxN Matrix with constraints coefficients
@param  b  Mx1 vector with constraints constants
@param  rho Rho value for the main algorithm
@param  IT  Maximum of iterations for the main algorithm
@return A Nx1 optimal solutions vector
*/

// P es NxN
// q es Nx1
// A es MxN
// b es Mx1

#ifdef _ADMM_SW_
void ADMM_wrapper(int IT, DataType (&c_ADMM)[M][P], DataType (&c_neg_ADMM)[M][P], DataType (&x)[N][P], DataType (&z)[M][P], DataType (&u)[M][P])
{
    ADMM(IT, q_neg_matlab, rho_neg_A_T_matlab, A_matlab, A_neg_matlab, c_ADMM, c_neg_ADMM, R_inv_matlab, x, z, u);
}
#else
#include "hls_stream.h"
void ADMM_wrapper(DataType (&c_ADMM)[M][P], DataType (&c_neg_ADMM)[M][P], DataType (&x)[N][P], DataType (&z)[M][P], DataType (&u)[M][P])
{
    ADMM<IT_hw>(q_neg_matlab, rho_neg_A_T_matlab, A_matlab, A_neg_matlab, c_ADMM, c_neg_ADMM, R_inv_matlab, x, z, u);
}

extern "C" {
void ADMM_wrapper_accel(hls::stream<axis_t> &in, hls::stream<axis_t> &out) {
#pragma HLS INTERFACE s_axilite port = return bundle = control
#pragma HLS INTERFACE axis port = in
#pragma HLS INTERFACE axis port = out

	static DataType x[N][P] = {{}};
    static DataType z[M][P] = {{}};
    static DataType u[M][P] = {{}};

	DataType c_ADMM_DMA[M_16][P];

#pragma HLS ARRAY_PARTITION variable = c_ADMM_DMA factor = 16 dim = 1 cyclic

	int j_limit = 512 / DataTypeSize;
	int i_limit = M_16 / j_limit;
	int i_limit_out = N / j_limit;
	converter_t converter;

load_c:
	for (int i = 0; i < i_limit; i++) {
	axis_t temp = in.read();
	for (int j = 0; j < j_limit; j++) {
			int high = j * DataTypeSize + DataTypeSize - 1;
			int low = j * DataTypeSize;
			int index = i * 16 + j;

			converter.i = temp.data.range(high, low);
			c_ADMM_DMA[index][0] = converter.d;
		}
	}
	DataType c_ADMM[M][P];
	DataType c_ADMM_neg[M][P];
	for(int i = 0; i < M; i++){
		c_ADMM[i][0]     =  c_ADMM_DMA[i][0];
		c_ADMM_neg[i][0] = -c_ADMM_DMA[i][0];
	}
	ADMM_wrapper(c_ADMM, c_ADMM_neg, x, z, u);

writex:
	for (int i = 0; i < i_limit_out; i++) {
		axis_t temp;
		for (int j = 0; j < j_limit; j++) {
			int high = j * DataTypeSize + DataTypeSize - 1;
			int low = j * DataTypeSize;

			converter.d = x[i * 16 + j][0];
			temp.data.range(high, low) = converter.i;
		}
		ap_uint<1> last = 0;
		if (i == i_limit_out - 1) {
		  last = 1;
		}
		temp.last = last;
		temp.keep = -1; // enabling all bytes
		out.write(temp);
	}

}
}
#endif
