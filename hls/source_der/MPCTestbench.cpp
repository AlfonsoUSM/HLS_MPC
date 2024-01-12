#include <iostream>
#include <fstream>
#include "system.hpp"
#include "mpc.hpp"
#include <iomanip>
#include <math.h>

using namespace std;
#define DISPLAY
#ifdef DENSE
char input_file_name[] = "MPC_der_dense_N2.bin";//"samples.bin";//
#else
char input_file_name[] = "MPC_der_sparse_N4.bin";
#endif
char output_file_name[]= "MPC_der_dense.csv";
double threshold = 1e-0;
// MatLab -> C++
// single -> float
// double -> double
typedef float sample_data_t;

int main(int argc, char *argv[]){

    cout << "Horizon : " << N_HOR << endl;
    // number of samples to read from samples file
    // large numbers can make the test bench really slow
    int nSamplestb = 2500;//2501;

//	if (argc!=2){
//	        cerr << "Must specify .bin\n";
//	        return 1;
//	}

//	ifstream samples(argv[1], ios::binary);
	ifstream samples(input_file_name, ios::binary);
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

	data_t A[N_SYS][N_SYS];
	data_t B[N_SYS][M_SYS];

    int aux = 0;
    sample_data_t aux2 = 0;

    // get number parameters from .bin file
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sN_SYS = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sM_SYS = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
	sP_SYS = aux;
    samples.read(reinterpret_cast<char*>(&aux), sizeof(uint8_t));
    sN_HOR = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint16_t));
	sN_QP = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint16_t));
	sM_QP = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint16_t));
	sIter = aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(uint16_t));
	sN = aux;
    aux = 0;

    // check if sN_HOR, sN_QP and sM_QP from .bin file match with N_HOR, N_QP and M_QP from system.hpp
    if (sN_HOR != N_HOR){
        cout << "N_HOR from sample file (" << sN_HOR << ") does not match N_HOR from system.hpp (";
        cout << N_HOR << ")" << endl;
        return 1;
    }
    if (sN_QP != N_QP){
        cerr << "N_QP from sample file (" << sN_QP << ")does not match N_QP from system.hpp (";
        cerr << N_QP << ")" << endl;
        return 1;
    }
    if (sM_QP != M_QP){
        cerr << "M_QP from sample file (" << sM_QP << ")does not match M_QP from system.hpp (";
        cerr << M_QP << ")" << endl;
        return 1;
    }

    // number of samples to read
	nSamplestb = (nSamplestb>sN)? sN : nSamplestb;

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

    data_t x[nSamplestb][N_SYS];
    data_t r[nSamplestb][P_SYS];
    data_t d[nSamplestb][D_SYS];
    data_t ref_u[nSamplestb][M_SYS];

    for (int sample=0; sample<nSamplestb; sample++){	// for each sample
        if (samples.peek() == EOF) break;	// exit if end of file
        for (int i=0; i<N_SYS; i++){		// load x
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            x[sample][i] = (data_t)aux2;
        }
        for (int i=0; i<P_SYS; i++){		// load r
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            r[sample][i] = (data_t)aux2;
        }
        for (int i=0; i<D_SYS; i++){		// load r
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            d[sample][i] = (data_t)aux2;
        }
        for (int i=0; i<M_SYS; i++){		// load expected u
        	samples.read(reinterpret_cast<char*>(&aux2), sizeof(sample_data_t));
            ref_u[sample][i] = aux2;
        }
    }

    ///////////////// Data from golden reference is ready ///////////////////

    cout << "Reference ready" << endl;

    data_t x0[N_SYS];
    data_t r0[P_SYS];
    data_t d0[D_SYS];
    data_t ref_u0[M_SYS];
    data_t u0[M_SYS] = {0};
//    data_t x1[N_SYS] = {0};

	data_t max_error = 0;
	data_t per_error;
	fstream fout;
	fout.open(output_file_name, ios::out);
	fout << "x0_1, x0_2, r0_1 r0_2 d0_1 d0_2 u0ref_1 u0ref_2, u0\n";

    for (int sample=0; sample<nSamplestb; sample++){	// for each sample
    	for (int i=0; i<N_SYS; i++){		// load x0
    	    x0[i] = x[sample][i];
    	    fout << x0[i] << ",";
    	}
    	for (int i=0; i<P_SYS; i++){		// load r0
    	    r0[i] = r[sample][i];
    	    fout << r0[i] << ",";
    	}
    	for (int i=0; i<D_SYS; i++){		// load d0
    	    d0[i] = d[sample][i];
    	    fout << d0[i] << ",";
    	}
        for (int i=0; i<M_SYS; i++){		// load ref u0
            ref_u0[i] = ref_u[sample][i];
        }
        mpc(x0, r0, d0, u0, 10);

#ifdef DISPLAY
    	cout << "sample number : " << sample << endl;
#endif
        for (int i=0; i<M_SYS; i++){		// error
        	double error = fabs((double)ref_u0[i] -(double)u0[i]);
        	double per_error = 100 * fabs((double)ref_u0[i] -(double)u0[i]) / fabs(ref_u0[i]);
        	fout << ref_u0[i] << "," << u0[i] << ",";
#ifdef DISPLAY
        	cout << "expected : " << ref_u0[i] << "\tresult : " << u0[i] << "\terror : " << error << "\t=>  " << per_error << " %" << endl;
#endif
        	if (fabs(ref_u0[i])>0.0002){
        		if (per_error > max_error){
        			max_error = per_error;
        		}
        	}
        }
        fout << "\n";
        cout << endl;
    }
    cout << "Max error : " << max_error << " %" << endl;

    samples.close();
	return 0;
}
