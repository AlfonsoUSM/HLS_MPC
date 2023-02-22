#include "linearSolversHW.h"


void minresHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zko[N_HOR], elem zk[N_HOR], elem iter, elem tol){
	// -----------------------------------------------------------------------

	int ip= 0;
	elem btemp = 0, beta, c=1, c_old=1, s=0, eta;
	elem norm_r, norm_r0, beta_old, rd, c_oold, alpha;
	elem s_oold, s_old, r1_hat, r1, r2, r3, ab_s;
	elem x[N_HOR], v[N_HOR], w[N_HOR], v_hat[N_HOR], v_old[N_HOR], w_old[N_HOR], w_oold[N_HOR], Av[N_HOR];
	// -------------------x = zko---------------------------------------------
	int i,k;
	initialize:
	for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
		x[k] = zko[k];
		v[k] = 0;
		w[k] = 0;
		w_old[k] = 0;
		Av[k] = 0;
	}
	// ----------------------Ak*x---------------------------------------------
	initializeAv:
	for (i = 0; i < N_HOR; ++i){
#pragma HLS UNROLL
		initializeAvN:
		for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
			Av[i] += Ak[i][k] * x[k];
		}
	}
	// ------------------r = b - Ak*x-----------------------------------------
	initializev_hat:
	for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
		v_hat[k] = bk[k]-Av[k] ;
	}
	// -----------------dw = r'*r---------------------------------------------
	btemp=0;
	initializebtemp:
	for (i = 0; i < N_HOR; ++i) {
#pragma HLS UNROLL
		btemp += v_hat[i] * v_hat[i];
	}
	// -------------tce = norm(r)---------------------------------------------
	beta = sqrt(btemp);
	eta=beta;
	norm_r = beta;
	norm_r0= beta;
    while_loop:while (ip<iter+1 && norm_r/norm_r0>tol){
#pragma HLS LOOP_TRIPCOUNT min=30 max=30 avg=30
		ip++;
		update1:
		for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
			v_old[k] = v[k];
			v[k] = v_hat[k]/beta;
		}
		// -------------------------Av=Ak*v------------------------------------
		//AvZeros:
		//for (k = 0; k < N_HOR; ++k){ ;}
		alpha = 0;
		Avcols:
		for (i = 0; i < N_HOR; ++i){
#pragma HLS UNROLL
			Av[i] = 0;
			Avrows:
			for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
				Av[i] += Ak[i][k] * v[k];
			}
			alpha += v[i] * Av[i];

		}
		// ---------------------dq=d'*q----------------------------------------
		rd = 0;
		v_hat:
		for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
			v_hat[k] = Av[k]-alpha*v[k]-beta*v_old[k];
			rd += v_hat[k] * v_hat[k];
		}
		beta_old=beta;
		// ----------------------------------------------------------------------

		beta = sqrt(rd);
		// ----------------------------------------------------------------------
		c_oold=c_old;
		c_old=c;                                            
        s_oold=s_old;                                        
        s_old=s; 
		r1_hat=c_old*alpha-c_oold*s_old*beta_old;           
        r1 =sqrt(r1_hat*r1_hat+beta*beta);
        r2 =s_old*alpha+c_oold*c_old*beta_old;              
        r3 =s_oold*beta_old;
		c=r1_hat/r1;                                        
        s=beta/r1;
		// ----------------------------------------------------------------------
        update2:
		for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
			w_oold[k] = w_old[k];
			w_old[k] = w[k];
			w[k]=(v[k]-r3*w_oold[k]-r2*w_old[k])/r1;
			x[k]=x[k]+c*eta*w[k];
		}
		if (s<0){
			ab_s=s*(-1);
		}else{
			ab_s=s;
		}
		norm_r = norm_r * ab_s;
		eta = -s*eta;
		// -----------r = b - Ak*x-----------------------------------------		
	}
	// -----------------------zk = x---------------------------------------
	for (k = 0; k < N_HOR; ++k){
		zk[k] = x[k];

	}
}


void cgradHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zko[N_HOR], elem zk[N_HOR], elem tol){
	// -------------------------------------------------------------------------
	int ip = 0; int k, i;
	elem dw, tce, dq, alpha, dl, beta, rd;
	elem x[N_HOR], x_old[N_HOR], r[N_HOR], d[N_HOR], q[N_HOR], rs[N_HOR], Ax[N_HOR];
	// -------------------x = zko-----------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		x[k] = zko[k];
		x_old[k] = zko[k];
		Ax[k] = 0;
	}
	// ----------------------Ak*x-----------------------------------------------
	for (i = 0; i < N_HOR; ++i){
		for (k = 0; k < N_HOR; ++k){
			Ax[i] += Ak[i][k] * x[k];
		}
	}
	// ------------------r = b - Ak*x-------------------------------------------
	dw=0;
	for (k = 0; k < N_HOR; ++k){
		r[k] = bk[k]-Ax[k] ;
		d[k] = r[k];
		dw += r[k] * r[k];
	}	
	// -------------tce = norm(r)-----------------------------------------------
	tce = sqrt(dw);
	while (ip<N_HOR && tce>tol){	
		// -------------------------q=Ak*d--------------------------------------
		for (k = 0; k < N_HOR; ++k){ q[k] = 0;}
		for (i = 0; i < N_HOR; ++i){
			for (k = 0; k < N_HOR; ++k){
				q[i] += Ak[i][k] * d[k];
			}
		}
		// ---------------------dq=d'*q-----------------------------------------
		dq = 0;
		for (i = 0; i < N_HOR; ++i){
			dq += d[i] * q[i];
		}
		alpha = dw/dq;
		// -----------r = b - Ak*x----------------------------------------------
		for (k = 0; k < N_HOR; ++k){
			x[k] = x[k]+alpha*d[k] ;
		}
		// ---------------------------------------------------------------------
		if (floor(ip/50)==ip/50){
			// ------------------Ak*x-------------------------------------------
			for (k = 0; k < N_HOR; ++k){ Ax[k] = 0;}
			for (i = 0; i < N_HOR; ++i){
				for (k = 0; k < N_HOR; ++k){
					Ax[i] += Ak[i][k] * x[k];
				}
			}
			// --------------------------r = b - Ak*x---------------------------
			for (k = 0; k < N_HOR; ++k){
				r[k] = bk[k]-Ax[k] ;
			}
		}else{
			// ----------------------r = r - alpha*q----------------------------
			for (k = 0; k < N_HOR; ++k){
				r[k] = r[k]-alpha*q[k] ;
			}
		}
		// ---------------------------------------------------------------------
		dl=dw;
		// -----------------------dw = r'*r-------------------------------------
		dw=0;
		for (i = 0; i < N_HOR; ++i) {
			dw += r[i] * r[i];
		}
		//----------------------------------------------------------------------
		beta=dw/dl;
		// -------------------d = r +beta*d-------------------------------------
		for (k = 0; k < N_HOR; ++k){
			d[k] = r[k]+beta*d[k] ;
		}
		// ---------------------tce=norm(x-zko)---------------------------------
		rd = 0;
		for (k = 0; k < N_HOR; ++k){
			rs[k] = x[k]-x_old[k] ;
			rd += rs[k] * rs[k];
		}
		// ---------------------------------------------------------------------
		tce = sqrt(rd);
		// -----------------------zko = x---------------------------------------
		for (k = 0; k < N_HOR; ++k){
			x_old[k] = x[k];
		}
		// ---------------------------------------------------------------------
		ip++;
	}
	// -----------------------zk = x--------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		zk[k] = x[k];
	}		
}


void cholHW(elem Ak[N_HOR][N_HOR], elem bk[N_HOR], elem zk[N_HOR]){
	// -------------------------------------------------------------------------

	int ik, tk, jv, iv, j, i, k, h, jz, t, ip,jp, kp;
	elem tol, sq, suma;
    elem x[N_HOR], y[N_HOR], Df[N_HOR], Rs[N_HOR], Tp[N_HOR], A[N_HOR*N_HOR], LT[N_HOR*N_HOR];
	

#if defined ELEM_DOUBLE
	elem eps = 2.220446049250313e-16;
#elif defined ELEM_FLOAT
	elem eps = 1.1920929e-07;
#else
#error Something is wrong with how elem is defined
#endif
	tol = N_HOR * eps;
	//--------Calcular L en Ak----------------------------------
	elem sqrtA = sqrt(Ak[0][0]);
	set_L0:for (int c=0; c<N_HOR; c++){
		elem aux = Ak[0][c]/sqrtA;
		if (Ak[0][0]<=tol){
			Ak[0][c] = 0;
		}
		else{
			Ak[0][c] = aux;
		}
	}

	set_L:for (int j=1; j<N_HOR; j++){
		for (int c=j; c<N_HOR; c++){
			for (int r=0; r<j; r++){
				Ak[j][c] -= Ak[r][j] * Ak[r][c];
			}
		}
		sqrtA = sqrt(Ak[j][j]);
		for (int c=j; c<N_HOR; c++){
		elem aux = Ak[j][c]/sqrtA;
			if (Ak[j][j]<=tol){
				Ak[j][c] = 0;
			}
			else{
				Ak[j][c] = aux;
			}
		}
	}

    set_A:for (int r = 0; r < N_HOR; ++r){
//#pragma HLS UNROLL
		for (int c = 0; c < N_HOR; ++c){
//#pragma HLS UNROLL
			if (r>c){
				A[N_HOR * r + c] = 0;
			}
			else{
				A[N_HOR * r + c] = Ak[r][c];
			}

		}
	}

	// ----------------TRANSPONER A --------------------------------------------

	Transp_A:for (i = 0; i < N_HOR; ++i){
#pragma HLS UNROLL
		for (jz = 0; jz < N_HOR; ++jz){
#pragma HLS UNROLL
			LT[N_HOR * i + jz]=A[N_HOR * jz + i];
		}			
	}
	//----------------CALCULO DE y----------------------------------------------
	suma=0;

	set_y:for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
		y[k] = 0;
	}
	for (i = 0; i < N_HOR; ++i){
#pragma HLS UNROLL
		for (jz = 0; jz < i; ++jz){
#pragma HLS UNROLL
			suma=suma+LT[N_HOR*i+jz]*y[jz];
		}	
		y[i]=(1/LT[N_HOR*i+i])*(bk[i]-suma);
		suma=0;
	}
	//----------------CALCULO DE x----------------------------------------------
	suma=0; tk=0;

	set_x:for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
		x[k] = 0;
	}
	for (i = 0; i < N_HOR; ++i){
#pragma HLS UNROLL
		tk=N_HOR-i-1;
		for (t = 0; t < i; ++t){
			suma=suma+A[N_HOR*tk-t+N_HOR-1]*x[N_HOR-t-1];
		}		
		x[N_HOR-i-1]=(1/LT[N_HOR*tk+tk])*(y[N_HOR-i-1]-suma);
		suma=0;
	}
	// -----------------------zk = x--------------------------------------------
	output:
	for (k = 0; k < N_HOR; ++k){
#pragma HLS UNROLL
		zk[k] = x[k];
	}	
}

