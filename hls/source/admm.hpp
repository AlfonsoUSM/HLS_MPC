#pragma once

#include "utils.hpp"

/*!
@brief  Alternating Method of Multipliers for solving quadratic convex optimization.
        Quadratic programing (QP) problem solution:
                    min 1/2 x'*Q*x + q'*x
                    s.t     A*x <= c
        Arrange the constraints into equality constraints for ADMM formulation, being finally:
                    min 1/2 x'*Q*x + q'*x + g(z)
                    s.t     A*x + z = c
                            z >= 0
        where g(z) is the indicator function of Z
            g(z) = 0    if z in Z
            g(z) = ∞    if any component of z not in Z

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

// Q es NxN
// q es Nx1
// A es MxN
// c es Mx1
#ifdef _ADMM_SW_
#include <stdlib.h>
template<int N, int M, int P, typename T>
void ADMM(int IT,           T (&q_neg)[N][P], T (&rho_neg_A_T)[N][M], 
          T (&A)[M][N],     T (&A_neg)[M][N], T (&c)[M][P],     
          T (&c_neg)[M][P], T (&R_inv)[N][N], 
          T (&x)[N][P],     T (&z)[M][P],     T (&u)[M][P]){
    // Optimization iterations
    for(int i = 0; i < IT; i++){
        x_min(q_neg, rho_neg_A_T, c_neg, R_inv, x, z, u);
        z_min(A_neg, c, x, z, u);
        u_update(A, c_neg, x, z, u);
    }
}
#else
template<int IT, int N, int M, int P, typename T>
void ADMM(T (&q_neg)[N][P], T (&rho_neg_A_T)[N][M], T (&A)[M][N],     T (&A_neg)[M][N], 
          T (&c)[M][P],     T (&c_neg)[M][P],       T (&R_inv)[N][N], 
          T (&x)[N][P],     T (&z)[M][P],           T (&u)[M][P]){
    // Optimization iterations
    for(int i = 0; i < IT; i++){
#pragma HLS PIPELINE
        x_min(q_neg, rho_neg_A_T, c_neg, R_inv, x, z, u);
        z_min(A_neg, c, x, z, u);
        u_update(A, c_neg, x, z, u);
    }
}
#endif

