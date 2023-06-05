#include <iostream>
#include <fstream>
#include "system.hpp"
#include "mpc.hpp"
#include <iomanip>
//#include <math.h>

using namespace std;
#define DISPLAY
char file_name[] = "MPC_motor_N16_double.bin"; //"MPC_motor_N4_double.bin";
double threshold = 1e-0;
// MatLab -> C++
// single -> float
// double -> double
typedef double sample_data_t;

int main(int argc, char *argv[]){

    // number of samples to read from samples file
    // large numbers can make the test bench really slow
    int nSamplestb = 500000;

//	if (argc!=2){
//	        cerr << "Must specify .bin\n";
//	        return 1;
//	}

//	ifstream samples(argv[1], ios::binary);
	ifstream samples(file_name, ios::binary);
	//check to see that the file was opened correctly:
	if (!samples.is_open()) {
	    cerr << "There was a problem opening the input file: ";
	    cerr << argv[1] << endl;
	    return 1;//exit or do additional error checking
	}

    int sN_HOR = 0;
    int sN_SYS = 0;
    int sM_SYS = 0;
    int sP_SYS = 0;
    int sN = 0;
    int sN_QP = 0;
    int sM_QP = 0;
    int sIter = 0;

    ///////////////// Read data from golden reference ///////////////////

	data_t xmin[N_SYS];
	data_t xmax[N_SYS];
	data_t umin[M_SYS];
	data_t umax[M_SYS];
//	data_t H_qp[N_QP][N_QP];
//	data_t h_qp[N_QP];
//	data_t C_qp[M_QP][N_QP];
//	data_t g[(2*N_QP)];
	data_t A[N_SYS][N_SYS];
	data_t B[N_SYS][M_SYS];

    int aux = 0;
    sample_data_t aux2 = 0;

    // get number parameters from .bin file
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sN_HOR = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sN_SYS = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sM_SYS = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sP_SYS = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint64_t));
	sN = aux;
	aux = 0;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sN_QP = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sM_QP = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint16_t));
	sIter = aux;
	aux = 0;

    // check if sN_HOR, sN_QP and sM_QP from .bin file match with N_HOR, N_QP and M_QP from system.hpp
    if (sN_HOR != N_HOR){
        cerr << "N_HOR from " << argv[1] << " (" << sN_HOR << ") does not match N_HOR from system.hpp (";
        cerr << N_HOR << ")" << endl;
        cout << "N_HOR from " << argv[1] << " (" << sN_HOR << ") does not match N_HOR from system.hpp (";
        cout << N_HOR << ")" << endl;
        return 1;
    }
    if (sN_QP != N_QP){
        cerr << "N_QP from " << argv[1] << " (" << sN_QP << ")does not match N_QP from system.hpp (";
        cerr << N_QP << ")" << endl;
        cout << "N_QP from " << argv[1] << " (" << sN_QP << ")does not match N_QP from system.hpp (";
        cout << N_QP << ")" << endl;
        return 1;
    }
    if (sM_QP != M_QP){
        cerr << "M_QP from " << argv[1] << " (" << sM_QP << ")does not match M_QP from system.hpp (";
        cerr << M_QP << ")" << endl;
        cout << "M_QP from " << argv[1] << " (" << sM_QP << ")does not match M_QP from system.hpp (";
        cout << M_QP << ")" << endl;
        return 1;
    }

    // number of samples to read
	nSamplestb = (nSamplestb>sN)? sN : nSamplestb;

    for (int r=0; r<N_SYS; r++){	// load xmin
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	xmin[r] = (data_t)aux2;
    }
    for (int r=0; r<N_SYS; r++){	// load xmax
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	xmax[r] = (data_t)aux2;
    }
    for (int r=0; r<M_SYS; r++){	// load umin
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	umin[r] = (data_t)aux2;
    }
    for (int r=0; r<M_SYS; r++){	// load umax
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	umax[r] = (data_t)aux2;
    }
    for (int r=0; r<N_QP; r++){	// load H_qp
    	for (int c=0; c<N_QP; c++){
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    		H_qp[r][c] = (data_t)aux2;
    	}
    }
    for (int r=0; r<N_QP; r++){	// load h_qp
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	h_qp[r] = (data_t)aux2;
    }
	for (int r=0; r<M_QP; r++){	// load C_qp
		for (int c=0; c<N_QP; c++){
	    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
			C_qp[r][c] = (data_t)aux2;
		}
    }
    for (int r=0; r<(2*N_QP); r++){	// load g
    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
    	g[r] = (data_t)aux2;
    }
	for (int r=0; r<N_SYS; r++){	// load A - Ignored in this testbench
		for (int c=0; c<N_SYS; c++){
	    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
			A[r][c] = (data_t)aux2;
		}
    }
	for (int r=0; r<N_SYS; r++){	// load B - Ignored in this testbench
		for (int c=0; c<M_SYS; c++){
	    	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
			B[r][c] = (data_t)aux2;
		}
    }

	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
	Rho = (data_t)aux2;

    for (int r=0; r<N_QP; r++){	// load R_inv
    	for (int c=0; c<N_QP; c++){
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
        	R_inv[r][c] = (data_t)aux2;
    	}
    }

    for (int r=0; r<N_QP; r++){	// load RhoMt_neg
    	for (int c=0; c<M_QP; c++){
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
        	RhoMt_neg[r][c] = (data_t)aux2;
    	}
    }

    data_t x[nSamplestb][N_SYS];
    data_t r[nSamplestb][P_SYS];
    data_t ref_u[nSamplestb][M_SYS];

    for (int sample=0; sample<nSamplestb; sample++){	// for each sample
        if (samples.peek() == EOF) break;	// exit if end of file
        for (int i=0; i<N_SYS; i++){		// load x
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            x[sample][i] = (data_t)aux2;
        }
//        for (int i=0; i<P_SYS; i++){		// load r
//        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
//            r[sample][i] = (data_t)aux2;
//        }
        for (int i=0; i<M_SYS; i++){		// load expected u
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            ref_u[sample][i] = aux2;
        }
    }

    ///////////////// Data from golden reference is ready ///////////////////

    data_t x0[N_SYS];
    data_t r0[P_SYS] = {0};
    data_t ref_u0[M_SYS];
    data_t u0[M_SYS] = {0};
    data_t x1[N_SYS] = {0};

	data_t theta[N_QP] = {0};
	data_t z_admm[M_QP] = {0};
	data_t u_admm[M_QP] = {0};

    for (int sample=0; sample<nSamplestb; sample++){	// for each sample
    	for (int i=0; i<N_SYS; i++){		// load x0
    	    x0[i] = x[sample][i];
    	}
        for (int i=0; i<M_SYS; i++){		// load ref u0
            ref_u0[i] = ref_u[sample][i];
        }
        mpc_sparse_admm_iteration(r0, x0, theta, z_admm, u_admm);
    	// estimate x1
    	for (int i=0; i<M_SYS ; i++){
    		u0[i] = theta[(N_HOR*N_SYS+N_SYS + i)];
    	}
    	// x1 =

        //cout << "sample number : " << sample << endl;
    	cout << "sample number : " << sample << endl;
        for (int i=0; i<M_SYS; i++){
        	double error = (double)ref_u0[i] -(double)u0[i];	// falta fabs()
        	cout << "expected : " << ref_u0[i] << "\tresult : " << u0[i] << "\terror : " << error << endl;

//            if ((error > threshold) | (error != error)){
//                errors++;
//#ifdef DISPLAY
//                cout << "sample number : " << sample << " " << i << endl;
//                cout << expected_u[i] << "\t";
//                cout << u0[i] << "\t";
//                cout << error;
//                cout << "\n";
//#endif
//            }
        }
        cout << endl;

    }

    samples.close();
	return 0;
}
