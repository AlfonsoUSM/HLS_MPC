#include "pdipSW.h"

void pdipSW(elem H[N_HOR][N_HOR], elem h[N_HOR], elem M[I_AUX][N_HOR], elem c[I_AUX], elem tk[N_HOR], int iterPDIP, int iterLs, elem tol){

	elem sgk = 0.5;
    elem bt = 0.99999;
    elem alp = 1;

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
    elem Dlk[I_AUX];
    elem Dlk_aux[I_AUX];
    elem Dsk[I_AUX];

    // MATLAB: tk=1*ones(n,1);
    for (int i=0; i<N_HOR; i++){
        tk[i] = 1;
    }
    // MATLAB: lk=0.5*ones(i,1); 
    // MATLAB: sk=0.5*ones(i,1);
    for (int i=0; i<I_AUX; i++){
        lk[i] = .5;
        sk[i] = .5;
    }

    // Precompute M'
	for (int r=0; r<N_HOR; r++){
		for(int c=0; c<I_AUX; c++){
			M_t[r][c] = M[c][r];
		}
	}
    // PDIP iterations
    for (int iter=0; iter<iterPDIP; iter++){
        // -------------- Build Ak -----------------------

        // MATLAB: RK=diag(lk./sk);
        // RK_diag is a vector with the diagonal
        for (int i=0; i<I_AUX; i++){
            RK_diag[i] = lk[i]/sk[i];
        }

        // MATLAB: Ak=H+M'*RK*M; 
        // Ak_aux = M'*RK
        for (int r=0; r<N_HOR; r++){
            for(int c=0; c<I_AUX; c++){
            	Ak_aux[r][c] = M_t[r][c]*RK_diag[c];
            }
        }
        // Ak = Ak_aux*M
        mmultSW (&Ak_aux[0][0], &M[0][0], &Ak[0][0], N_HOR, I_AUX, N_HOR);
        // Ak = Ak+H
        for (int r=0; r<N_HOR; r++){
            for(int c=0; c<N_HOR; c++){
            	Ak[r][c] = Ak[r][c]+H[r][c];
            }
        }

        // -------------- Build bk -----------------------

        // MATLAB: muk=(lk'*sk)/i;
        // muk = (lk'*sk)
        muk = 0;
        for (int i=0; i<I_AUX; i++){
        	muk += lk[i] * sk[i];
        }
        // muk = muk/i
        muk = muk/I_AUX;

        // MATLAB: HK=-H*tk-h-M'*lk;
        // HK_aux1 = H*tk
        mmultSW (&H[0][0], &tk[0], &HK_aux1[0], N_HOR, N_HOR, 1);
        // HK_aux2 = M'*lk
        mmultSW (&M_t[0][0], &lk[0], &HK_aux2[0], N_HOR, I_AUX, 1);
        // HK = -HK_aux1 -h -HK_aux2
        for (int i=0; i<N_HOR; i++){
            HK[i] = -HK_aux1[i] - h[i] - HK_aux2[i];
        }

        // MATLAB: GK=-M*tk+c-sk;
        // GK_aux = M*tk
        // GK = -GK_aux +c -sk (in set_GK_bk_aux1_TK)
        mmultSW (&M[0][0], &tk[0], &GK_aux[0], I_AUX, N_HOR, 1);
        
        // MATLAB: TK=-lk.*sk+sgk*muk*em;
        // here TK is divided by lk so that it's not 
        // necessary to divide by it later
        // TK = -sk + sgk*muk/lk (in set_GK_bk_aux1_TK)

        // MATLAB: bk= HK+M'*RK*(GK-TK./lk); // here TK is TK./lk
        // bk_aux = (GK+TK)
        for (int i=0; i<I_AUX; i++){
        	GK[i] = -GK_aux[i] + c[i] - sk[i];
        	TK[i] = -sk[i]+sgk*muk/lk[i];
            bk_aux[i] = GK[i]-TK[i];
        }
        // remember that Ak_aux = M'*Rk
        // bk = Ak_aux*bk_aux
        mmultSW (&Ak_aux[0][0], &bk_aux[0], &bk[0], N_HOR, I_AUX, 1);
        // bk = HK+bk
        for (int i=0; i<N_HOR; i++){
            bk[i] += HK[i];
        }

        // --------- Solve the system Ak*zk=bk ----------

#if LS == MINRES 
minresSW(Ak, bk, zko, zk, iterLs, tol);
#elif LS == CGRAD
    cgradSW(Ak, bk, zko, zk, tol);
#elif LS == CHOL
    cholSW(&Ak[0][0], &bk[0], &zk[0]);
#else
#error Something is wrong with how the linear solver is defined in specs.h
#endif

        // ----- Calculate Delta_lk , Delta_sk and find optimal alp ------
        //                              and
        // ----------- -----------Find best alp --------------------------
        
        // MATLAB: Dlk=-RK*(GK-TK./lk-M*zk); // here TK is TK./lk
        // Dlk_aux = M*zk
        mmultSW (&M[0][0], &zk[0], &Dlk_aux[0], I_AUX, N_HOR, 1);
        
        // The search for the best alp is split in two,
        // the first half is done while making Dlk and 
        // the second one is done while making Dsk.

        alp = 1;
        // Dlk = RK_diag * (GK-TK-Dlk_aux1)
        for (int i=0; i<I_AUX; i++){
            Dlk[i] = -RK_diag[i] * (GK[i]-TK[i] - Dlk_aux[i]);
            elem Dlk_i = (Dlk[i]<0) ? Dlk[i] : 0;
            
            // MATLAB: alp_lk = max(alp_lk,max(-Dlk(idx)./lk(idx)));
            Dlk_i = -Dlk_i/lk[i];
            alp = (alp>Dlk_i)? alp : Dlk_i;
        }

        // MATLAB: Dsk=TK./lk-RKI*Dlk; // here TK is TK./lk
        // RKI is 1/RK_diag
        for (int i=0; i<I_AUX; i++){
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
        for (int i=0; i<N_HOR; i++){
            tk[i] = tk[i] +alp*zk[i];
            zko[i] = zk[i];
        }

        // MATLAB: lk=lk+alp*Dlk;
        // MATLAB: sk=sk+alp*Dsk;
        for (int i=0; i<I_AUX; i++){
            lk[i] = lk[i] +alp*Dlk[i];
            sk[i] = sk[i] +alp*Dsk[i];
        }
    }

}



