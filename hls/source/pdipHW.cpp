#include "pdipHW.h"

void pdipHW(elem H[N_HOR][N_HOR], elem h[N_HOR], elem M[I_AUX][N_HOR], elem c[I_AUX], elem tk[N_HOR], int iterPDIP, int iterLs, elem tol){

	elem sgk = 0.5;
    elem bt = 0.99999;
    elem alp = 1;

    elem Wk[I_AUX][I_AUX];
    elem RK_diag[I_AUX];
	elem Ak[N_HOR][N_HOR];
	elem Ak_aux[N_HOR][I_AUX];
	elem M_t[N_HOR][I_AUX];
	elem bk_aux[I_AUX];
	elem HK_aux1[N_HOR], HK_aux2[N_HOR];
	elem zk[N_HOR];
    elem zko[N_HOR] = {0}; 


    elem bk[N_HOR];
    elem muk;
    elem muk_aux[I_AUX];
    elem HK[N_HOR], GK[I_AUX];
    elem GK_aux[I_AUX];
    elem lk[I_AUX], sk[I_AUX];
    elem TK[I_AUX];
    elem TK_aux;
    elem Dlk[I_AUX];
    elem Dlk_aux[I_AUX];
    elem Dsk[I_AUX];

#pragma HLS ARRAY_PARTITION variable=Ak_aux complete dim=2
#pragma HLS ARRAY_PARTITION variable=M complete dim=1

#pragma HLS ARRAY_PARTITION variable=H complete dim=2
#pragma HLS ARRAY_PARTITION variable=tk complete dim=1

#pragma HLS ARRAY_PARTITION variable=M_t complete dim=2
#pragma HLS ARRAY_PARTITION variable=lk complete dim=1

#pragma HLS ARRAY_PARTITION variable=Ak complete dim=2

#pragma HLS ARRAY_PARTITION variable=bk_aux1 complete dim=1

#pragma HLS ARRAY_PARTITION variable=RK_diag complete dim=1

#pragma HLS ARRAY_PARTITION variable=sk complete dim=1

#pragma HLS ARRAY_PARTITION variable=zk complete dim=1

//#pragma HLS ARRAY_PARTITION variable=Dlk complete dim=1
//#pragma HLS ARRAY_PARTITION variable=Dsk complete dim=1

    // MATLAB: tk=1*ones(n,1);
    set_tk:
    for (int i=0; i<N_HOR; i++){
#pragma HLS UNROLL
        tk[i] = 1;
    }
    // MATLAB: lk=0.5*ones(i,1); 
    // MATLAB: sk=0.5*ones(i,1);
    set_lksk:
    for (int i=0; i<I_AUX; i++){
#pragma HLS UNROLL
        lk[i] = .5;
        sk[i] = .5;
    }

    // Precompute M'
	transpose_M:
	for (int r=0; r<N_HOR; r++){
		#pragma HLS pipeline
		for(int c=0; c<I_AUX; c++){
			M_t[r][c] = M[c][r];
		}
	}
    // PDIP iterations
    mainForLoop:
    for (int iter=0; iter<iterPDIP; iter++){
#pragma HLS LOOP_TRIPCOUNT min=20 max=20 avg=20

        // -------------- Build Ak -----------------------

        // MATLAB: RK=diag(lk./sk);
        // RK_diag is a vector with the diagonal
    	set_RK_diag:
        for (int i=0; i<I_AUX; i++){
			#pragma HLS unroll
            RK_diag[i] = lk[i]/sk[i];
        }

        // MATLAB: Ak=H+M'*RK*M; 
        // Ak_aux = M'*RK
        set_Ak_aux:
        for (int r=0; r<N_HOR; r++){
			#pragma HLS pipeline
            for(int c=0; c<I_AUX; c++){
            	Ak_aux[r][c] = M_t[r][c]*RK_diag[c];
            }
        }
        // Ak = Ak_aux*M
        mmult_N_I_N (Ak_aux, M, Ak);
        // Ak = Ak+H
        add_Ak_H:
        for (int r=0; r<N_HOR; r++){
			#pragma HLS pipeline
            for(int c=0; c<N_HOR; c++){
            	Ak[r][c] = Ak[r][c]+H[r][c];
            }
        }

        // -------------- Build bk -----------------------

        // MATLAB: muk=(lk'*sk)/i;
        // muk = (lk'*sk)
        muk = 0;
        set_muk_aux:
        for (int i=0; i<I_AUX; i++){
			#pragma HLS unroll
        	muk_aux[i] = lk[i] * sk[i];
        	muk += muk_aux[i];
        }
        // muk = muk/i
        muk = muk/I_AUX;

        // MATLAB: HK=-H*tk-h-M'*lk;
        // HK_aux1 = H*tk
        mmult_N_N_1 (H, tk, HK_aux1);
        // HK_aux2 = M'*lk
        mmult_N_I_1 (M_t, lk, HK_aux2);
        // HK = -HK_aux1 -h -HK_aux2
        set_HK:
        for (int i=0; i<N_HOR; i++){
			#pragma HLS pipeline
            HK[i] = -HK_aux1[i] - h[i] - HK_aux2[i];
        }

        // MATLAB: GK=-M*tk+c-sk;
        // GK_aux = M*tk
        // GK = -GK_aux +c -sk (in set_GK_bk_aux1_TK)
        mmult_I_N_1 (M, tk, GK_aux);
        
        // MATLAB: TK=-lk.*sk+sgk*muk*em;
        // here TK is divided by lk so that it's not 
        // necessary to divide by it later
        // TK = -sk + sgk*muk/lk (in set_GK_bk_aux1_TK)

        // MATLAB: bk= HK+M'*RK*(GK-TK./lk); // here TK is TK./lk
        // bk_aux = (GK+TK)
        set_GK_bk_aux1_TK:
        for (int i=0; i<I_AUX; i++){
			#pragma HLS pipeline
        	GK[i] = -GK_aux[i] + c[i] - sk[i];
        	TK[i] = -sk[i]+sgk*muk/lk[i];
            bk_aux[i] = GK[i]-TK[i];
        }
        // remember that Ak_aux = M'*Rk
        // bk = Ak_aux*bk_aux
        mmult_N_I_1 (Ak_aux, bk_aux, bk);
        // bk = HK+bk
        set_bk:
        for (int i=0; i<N_HOR; i++){
			#pragma HLS pipeline
            bk[i] += HK[i];
        }

        // --------- Solve the system Ak*zk=bk ----------

#if defined CHOL
        cholHW(Ak, bk, zk);
#elif defined MINRES
        minresHW(Ak, bk, zko, zk, iterLs, tol);
#elif defined CGRAD
        cgradHW(Ak, bk, zko, zk, tol);
#else
#error Something is wrong with how the linear solver is defined in specs.h
#endif

        // ----- Calculate Delta_lk , Delta_sk and find optimal alp ------
        //                              and
        // ----------- -----------Find best alp --------------------------
        
        // MATLAB: Dlk=-RK*(GK-TK./lk-M*zk); // here TK is TK./lk
        // Dlk_aux = M*zk
        mmult_I_N_1(M, zk, Dlk_aux);
        
        // The search for the best alp is split in two,
        // the first half is done while making Dlk and 
        // the second one is done while making Dsk.

        alp = 1;
        // Dlk = RK_diag * (GK-TK-Dlk_aux1)
        set_Dlk:
        for (int i=0; i<I_AUX; i++){
			#pragma HLS pipeline
            Dlk[i] = -RK_diag[i] * (GK[i]-TK[i] - Dlk_aux[i]);
            elem Dlk_i = (Dlk[i]<0) ? Dlk[i] : 0;
            
            // MATLAB: alp_lk = max(alp_lk,max(-Dlk(idx)./lk(idx)));
            Dlk_i = -Dlk_i/lk[i];
            alp = (alp>Dlk_i)? alp : Dlk_i;
        }

        // MATLAB: Dsk=TK./lk-RKI*Dlk; // here TK is TK./lk
        // RKI is 1/RK_diag
        set_Dsk:
        for (int i=0; i<I_AUX; i++){
			#pragma HLS pipeline
        	Dsk[i] = TK[i] - Dlk[i]/RK_diag[i];
        	elem Dsk_i = (Dsk[i]<0) ? Dsk[i] : 0;
        	
            // MATLAB: alp_sk = max(alp_sk,max(-Dsk(idx)./sk(idx)));
            Dsk_i = -Dsk_i/sk[i];
        	alp = (alp>Dsk_i)? alp : Dsk_i;
        }

        // MATLAB: alp=bt/max(alp_lk,alp_sk);
        // Here alp is already max(alp_lk,alp_sk)
        alp = bt/alp;

        // ----------- Prepare the next iteration ---------

        // MATLAB: tk=tk+alp*zk;
        // MATLAB: zko=zk;
        update_tk_zko:
        for (int i=0; i<N_HOR; i++){
			#pragma HLS PIPELINE
            tk[i] = tk[i] +alp*zk[i];
            zko[i] = zk[i];
        }

        // MATLAB: lk=lk+alp*Dlk;
        // MATLAB: sk=sk+alp*Dsk;
        update_lk_sk:
        for (int i=0; i<I_AUX; i++){
			//#pragma HLS pipeline
			#pragma HLS unroll
            lk[i] = lk[i] +alp*Dlk[i];
            sk[i] = sk[i] +alp*Dsk[i];
        }
    }

}



