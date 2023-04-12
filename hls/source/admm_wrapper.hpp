#pragma once

#include "ADMM.hpp"
#include "init_der.hpp"
#define N 16
#define M 74
#define P 1
#define M_16 80
#define N2 256 // N*N

#ifdef _ADMM_SW_
typedef float DataType;
extern void ADMM_wrapper(int IT, DataType (&c_ADMM)[M][P], DataType (&c_neg_ADMM)[M][P], DataType (&x)[N][P], DataType (&z)[M][P], DataType (&u)[M][P]);

#else
#include "ap_axi_sdata.h"
#include "ap_int.h"
#include "ap_fixed.h"
#include <inttypes.h>

#define IT_hw 1
#define rho_hw 62.963413f
#define DWIDTH 512

typedef float DataType;
//typedef ap_fixed<16, 4> DataType;
typedef ap_axiu<DWIDTH, 0, 0, 0> axis_t;

typedef ap_uint<512> uint512_t;

const int DataTypeSize = sizeof(DataType) * 8;

typedef ap_uint<DataTypeSize> DataTypeInt;

typedef union converter {
  DataType d;
  uint32_t i;
} converter_t;

extern void ADMM_wrapper(DataType (&c_ADMM)[M][P], DataType (&c_neg_ADMM)[M][P], DataType (&x)[N][P], DataType (&z)[M][P], DataType (&u)[M][P]);
#endif
