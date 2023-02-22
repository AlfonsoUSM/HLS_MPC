#include "linearSolversSW.h"

void minresSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem iter, elem tol){	
	// -------------------------------------------------------------------------
	int ip= 0;
	elem btemp = 0, beta, c=1, c_old=1, s=0, eta;
	elem norm_r, norm_r0, beta_old, rd, c_oold, alpha;
	elem s_oold, s_old, r1_hat, r1, r2, r3, ab_s;
	elem x[N_HOR], v[N_HOR], w[N_HOR], v_hat[N_HOR], v_old[N_HOR], w_old[N_HOR], w_oold[N_HOR], Av[N_HOR];
	// -------------------x = zko-----------------------------------------------
	int i,k;
	for (k = 0; k < N_HOR; ++k){
		x[k] = zko[k];
		v[k] = 0;
		w[k] = 0;
		w_old[k] = 0;
		Av[k] = 0;
	}
	// ----------------------Ak*x-----------------------------------------------
	for (i = 0; i < N_HOR; ++i){
		for (k = 0; k < N_HOR; ++k){
			Av[i] += Ak[N_HOR * i + k] * x[k];
		}
	}
	// ------------------r = b - Ak*x-------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		v_hat[k] = bk[k]-Av[k] ;
	}
	// -----------------dw = r'*r-----------------------------------------------
	btemp=0;
	for (i = 0; i < N_HOR; ++i) {
		btemp += v_hat[i] * v_hat[i];
	}
	// -------------tce = norm(r)-----------------------------------------------
	beta = sqrt(btemp);
	eta=beta;
	norm_r = beta;
	norm_r0= beta;
	while (ip<iter+1 && norm_r/norm_r0>tol){
		ip++;
		for (k = 0; k < N_HOR; ++k){
			v_old[k] = v[k];
			v[k] = v_hat[k]/beta;
		}
		// -------------------------Av=Ak*v-------------------------------------
		for (k = 0; k < N_HOR; ++k){ Av[k] = 0;}
		for (i = 0; i < N_HOR; ++i){
			for (k = 0; k < N_HOR; ++k){
				Av[i] += Ak[N_HOR * i + k] * v[k];
			}
		}
		// ---------------------dq=d'*q-----------------------------------------
		alpha = 0;
		for (i = 0; i < N_HOR; ++i){
			alpha += v[i] * Av[i];
		}
		for (k = 0; k < N_HOR; ++k){
			v_hat[k] = Av[k]-alpha*v[k]-beta*v_old[k];
		}
		beta_old=beta;
		// ---------------------------------------------------------------------
		rd = 0;
		for (i = 0; i < N_HOR; ++i) {
			rd += v_hat[i] * v_hat[i];
		}
		beta = sqrt(rd);
		// ---------------------------------------------------------------------
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
		// ---------------------------------------------------------------------
		for (k = 0; k < N_HOR; ++k){
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
		// -----------r = b - Ak*x----------------------------------------------
	}
	// -----------------------zk = x--------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		zk[k] = x[k];
	}
}

void cgradSW(elem* Ak, elem* bk, elem* zko, elem* zk, elem tol){	
	// -------------------------------------------------------------------------
	int ip = 0; int k, i;
	elem dw, tce, dq, alpha, dl, beta, rd;
	elem x[N_HOR], r[N_HOR], d[N_HOR], q[N_HOR], rs[N_HOR], Ax[N_HOR];
	// -------------------x = zko-----------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		x[k] = zko[k];
		Ax[k] = 0;
	}
	// ----------------------Ak*x-----------------------------------------------
	for (i = 0; i < N_HOR; ++i){
		for (k = 0; k < N_HOR; ++k){
			Ax[i] += Ak[N_HOR * i + k] * x[k];
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
				q[i] += Ak[N_HOR * i + k] * d[k];
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
					Ax[i] += Ak[N_HOR * i + k] * x[k];
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
			rs[k] = x[k]-zko[k] ;
			rd += rs[k] * rs[k];
		}
		// ---------------------------------------------------------------------
		tce = sqrt(rd);
		// -----------------------zko = x---------------------------------------
		for (k = 0; k < N_HOR; ++k){
			zko[k] = x[k];
		}
		// ---------------------------------------------------------------------
		ip++;
	}
	// -----------------------zk = x--------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		zk[k] = x[k];
	}		
}

void cholSW(elem* Ak, elem* bk, elem* zk){

	// -------------------------------------------------------------------------
	int ms = ceil((N_HOR/2))*ceil((N_HOR/2)), ik, tk, iv, j, i, k, h, jz, t, ip,jp, kp;
	elem tol, sq, suma;
	elem x[N_HOR], y[N_HOR], Df[N_HOR], Rs[N_HOR], Tp[N_HOR], Ts[ms], A[N_HOR*N_HOR], LT[N_HOR*N_HOR];
	
	// ----------------------B=Ak-----------------------------------------------
    int jv;
    for (iv = 0; iv < N_HOR; ++iv){
    	for (jv = 0; jv < N_HOR; ++jv){
    		A[N_HOR * iv + jv] = Ak[N_HOR * iv + jv] ;
    	}
    }
	
	// ------------------Obtener triangular superior de A-----------------------
	//for (iv = 1; iv < N_HOR; ++iv){
	//	for (jv = 0; jv < iv; ++jv){
	//		A[N_HOR * iv + jv]=0;
	//	}			
	//}
#if defined ELEM_DOUBLE
	elem eps = 2.220446049250313e-16;
#elif defined ELEM_FLOAT
	elem eps = 1.1920929e-07;
#else
#error Something is wrong with how elem is defined
#endif
	tol = N_HOR * eps;
	//--------------------------------------------------------------------------
	if (A[0] <= tol){
		for (iv = 0; iv < N_HOR; ++iv){
			A[iv] = 0;
		}
	}else{
		sq=sqrt(A[0]);
		for (iv = 0; iv < N_HOR; ++iv){
			A[iv] = A[iv]/sq;
		}
	}
	// --------------Constrir Tp y Ts-------------------------------------------
	j=0;
	for (j=1; j<N_HOR; ++j){
		//---------------------Tp-----------------------------------------------	 
		for (k=0; k<j; ++k){
			Tp[k]=A[k*N_HOR+j];
		}
		//---------------------Ts-----------------------------------------------
		ik=0;
		for (i=0; i<j; ++i){
			for (k=0; k<N_HOR-j; ++k){
				Ts[ik]=A[i*N_HOR + k + j];	
				ik=ik+1;
			}
		}
		//---------------------Rs-----------------------------------------------
		for (k = 0; k < N_HOR; ++k){ Rs[k] = 0;}		
		//for (ip = 0; ip < 1; ++ip){
			ip=0;
			for (jp = 0; jp < N_HOR-j; ++jp){
				for (kp = 0; kp < j; ++kp){
					Rs[(N_HOR-j) * ip + jp] += Tp[j * ip + kp] * Ts[(N_HOR-j) * kp + jp];
				}
			}
		//}
		//---------------------Df-----------------------------------------------
		for (h = 0; h < N_HOR-j; ++h){
			Df[h] = A[N_HOR*j + h+j];
		}
		//----------------------------------------------------------------------					
		for (h = 0; h < N_HOR-j; ++h){
			A[N_HOR*j + h+j] = Df[h]- Rs[h];
		}
		//----------------------------------------------------------------------
		if (A[N_HOR*j+j] <= tol){
			for (i = j; i < N_HOR; ++i){
				A[i] = 0;
			}
		}else{
			sq=sqrt(A[N_HOR*j+j]);
			for (i = 0; i < N_HOR-j; ++i){
				A[N_HOR*j+j+i] = A[N_HOR*j+j+i]/sq;
			}
		}		
	}	
	// ----------------TRANSPONER A --------------------------------------------
	for (i = 0; i < N_HOR; ++i){
		for (jz = 0; jz < N_HOR; ++jz){
			LT[N_HOR * i + jz]=A[N_HOR * jz + i];
		}			
	}
	//----------------CALCULO DE y----------------------------------------------
	suma=0;
	for (k = 0; k < N_HOR; ++k){ y[k] = 0;}
	for (i = 0; i < N_HOR; ++i){
		for (jz = 0; jz < i; ++jz){
			suma=suma+LT[N_HOR*i+jz]*y[jz];
		}	
		y[i]=(1/LT[N_HOR*i+i])*(bk[i]-suma);
		suma=0;
	}
	//----------------CALCULO DE x----------------------------------------------
	suma=0; tk=0;
	for (k = 0; k < N_HOR; ++k){ x[k] = 0;}
	for (i = 0; i < N_HOR; ++i){
		tk=N_HOR-i-1;
		for (t = 0; t < i; ++t){
			suma=suma+A[N_HOR*tk-t+N_HOR-1]*x[N_HOR-t-1];
		}		
		x[N_HOR-i-1]=(1/LT[N_HOR*tk+tk])*(y[N_HOR-i-1]-suma);
		suma=0;
	}
	// -----------------------zk = x--------------------------------------------
	for (k = 0; k < N_HOR; ++k){
		zk[k] = x[k];
    }


}