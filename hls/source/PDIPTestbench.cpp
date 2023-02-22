#include <iostream>
#include <fstream>
#include "specs.h"
#include "pdipHW.h"
#include <iomanip>

using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "PDIP testbench" << endl;
    elem H[N_HOR][N_HOR];
    elem h_tilde[N_HOR];
    elem M_hat[I_AUX][N_HOR];
    elem c_hat[I_AUX];
    double expectedResult[N_HOR];
    elem u_tilde_star[N_HOR];

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
    double threshold = 1e-4;
    // tolerance for the linear solver in PDIP algorithm
    elem tol = 1e-9;
    // PDIP iterations
    int iterPDIP = 20;
    // MINRES iterations
    int iterLS = 30;

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


    for (int sample=0; sample<nSamplestb; sample++){
    	// exit if end of file
        if (samples.peek() == EOF) break;
        // load h_tilde
        for (int i=0; i<N_HOR; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
        	h_tilde[i] = (elem)aux;
        }
        // load c_hat
        for (int i=0; i<I_AUX; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
        	c_hat[i] = (elem)aux;
        }
        // load expectedResult
        for (int i=0; i<N_HOR; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            expectedResult[i] = aux;
        }

        pdipHW(H, h_tilde, M_hat, c_hat, u_tilde_star, iterPDIP, iterLS, tol);


        //cout << "sample number : " << sample << endl;
        for (int i=0; i<N_HOR; i++){
        	double error = fabs(expectedResult[i] -(double)u_tilde_star[i]);
        	//cout << error<< "\n";
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << u_tilde_star[i] << "\t";
                cout << error;
                cout << "\n";
#endif

            }
        }
        
    }

    samples.close();

    // print important data
    cout << "Finished " << nSamplestb << " QP problems" << endl;
    cout << "Number of differences between expected and calculated " << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;

    return errors;
}
