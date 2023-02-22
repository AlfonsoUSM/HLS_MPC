#include <iostream>
#include <fstream>
#include "specs.h"
#include "mpcHW.h"
#include <iomanip>

using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "MPC testbench" << endl;
	elem x[N_SYS];
	elem r[P_SYS];
	elem xmin[N_SYS];
	elem xmax[N_SYS];
	elem umin[M_SYS];
	elem umax[M_SYS];
	elem Acal[N_SYS*N_HOR][N_SYS];
	elem AcalQOcal[N_SYS][N_HOR];
	elem H[N_HOR][N_HOR];
	elem M_hat[N_SYS*N_HOR*2+2*N_HOR][N_HOR];
	elem L_invLast[(N_SYS+P_SYS)];
    elem A[N_SYS][N_SYS];
    elem B[N_SYS][M_SYS];
	double expectedResult[P_SYS];
	elem u[M_SYS];
	int iterPDIP = 20;
	int iterLS = 30;
	int errors = 0;

    if (argc!=2){
        cerr << "Must specify .bin\n";
        return 1;
    }

    ifstream samples(argv[1], ios::binary);
    //check to see that the file was opened correctly:
    if (!samples.is_open()) {
        cerr << "There was a problem opening the input file: ";
        cerr << argv[1] << endl;
        return 1;//exit or do additional error checking
    }

    // threshold for the difference between the calculated results and the correct ones.
    double threshold = 1e-0;
    // tolerance for the linear solver in PDIP algorithm
    elem tol = 1e-9;

    int nSamples = 0;
    int samplesI = 0;
    int samplesN = 0;

    // number of samples to read from samples file
    // large numbers can make the test bench really slow
    int nSamplestb = 500000;

    // aux is used to temporary store the value before cast
    double aux;

    // get number N, I and number of samples
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	samplesN = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	samplesI = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	nSamples = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	iterPDIP = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	iterLS = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	tol = (int)aux;

    // check if N and I match with N_HOR and I_AUX from specs.h
    if (N_HOR != samplesN){
        cerr << "N_HOR from " << argv[1] << " (" << samplesN << ") does not match N_HOR from specs.h (";
        cerr << N_HOR << ")" << endl;
        return 1;
    }
    if (I_AUX != samplesI){
        cerr << "I_AUX from " << argv[1] << " (" << samplesI << ")does not match I_AUX from specs.h (";
        cerr << I_AUX << ")" << endl;
        return 1;
    }

    nSamplestb = (nSamplestb>nSamples)? nSamples : nSamplestb;

    // load xmin
    for (int r=0; r<N_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	xmin[r] = (elem)aux;
    }
    // load xmax
    for (int r=0; r<N_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	xmax[r] = (elem)aux;
    }
    // load umin
    for (int r=0; r<M_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	umin[r] = (elem)aux;
    }
    // load umax
    for (int r=0; r<M_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	umax[r] = (elem)aux;
    }
    // load Acal
    for (int r=0; r<N_SYS*N_HOR; r++){
    	for (int c=0; c<N_SYS; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		Acal[r][c] = (elem)aux;
    	}
    }
    // load AcalQOcal
    for (int r=0; r<N_SYS; r++){
    	for (int c=0; c<N_HOR; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		AcalQOcal[r][c] = (elem)aux;
    	}
    }
    // load H
    for (int r=0; r<N_HOR; r++){
    	for (int c=0; c<N_HOR; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		H[r][c] = (elem)aux;
    	}
    }
    // load M_hat
	for (int r=0; r<I_AUX; r++){
		for (int c=0; c<N_HOR; c++){
			samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			M_hat[r][c] = (elem)aux;
		}
    }
    // load L_last
	for (int r=0; r<N_SYS+P_SYS; r++){
		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
		L_invLast[r] = (elem)aux;
    }
    // load A - Ignored in this testbench
	for (int r=0; r<N_SYS; r++){
		for (int c=0; c<N_SYS; c++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			A[r][c] = (elem)aux;
		}
    }
    // load B - Ignored in this testbench
	for (int r=0; r<N_SYS; r++){
		for (int c=0; c<M_SYS; c++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			B[r][c] = (elem)aux;
		}
    }


    for (int sample=0; sample<nSamplestb; sample++){
    	// exit if end of file
        if (samples.peek() == EOF) break;
        // load x
        for (int i=0; i<N_SYS; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            x[i] = (elem)aux;
        }
        // load r
        for (int i=0; i<P_SYS; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            r[i] = (elem)aux;
        }
        // load expectedResult
        for (int i=0; i<P_SYS; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            expectedResult[i] = aux;
        }

        mpcHW(x, r, xmin, xmax, umin, umax, Acal, AcalQOcal, H, M_hat, L_invLast, iterPDIP, iterLS, u);

        //cout << "sample number : " << sample << endl;
        for (int i=0; i<M_SYS; i++){
        	double error = fabs(expectedResult[i] -(double)u[i]);
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << u[i] << "\t";
                cout << error;
                cout << "\n";
#endif

            }
        }

    }

    samples.close();

    // print important data
    cout << "Finished processing " << nSamplestb << " samples" << endl;
    cout << "Number of differences between expected and calculated " << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;

    return errors;
}
