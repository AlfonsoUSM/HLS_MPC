#include <iostream>
#include <fstream>
#include "specs.h"
#include "linearSolversHW.h"
#include <iomanip>



using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "Linear Solver testbench" << endl;
	elem A[N_HOR][N_HOR];
	elem b[N_HOR];
	elem z[N_HOR];
	double expectedResult[N_HOR];

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
    double threshold = 1e-4;
    // tolerance for the linear solvers MINRES and CGRAD
    elem tol = 1e-9;

    int nSamples = 0;
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
    nSamples = (int)aux;
    samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    iterLS = (int)aux;
	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
	tol = (int)aux;

    nSamplestb = (nSamplestb>nSamples)? nSamples : nSamplestb;

    // check if N match with N_HOR and I_QP from specs.h
    if (N_HOR != samplesN){
        cerr << "N_HOR from " << argv[1] << " (" << samplesN << ") does not match N_HOR from specs.h (";
        cerr << N_HOR << ")" << endl;
        return 1;
    }

    for (int sample=0; sample<nSamplestb; sample++){
    	// exit if end of file
        if (samples.peek() == EOF) break;
        // load A
        for (int r=0; r<N_HOR; r++){
        	for (int c=0; c<N_HOR; c++){
                samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
                A[r][c] = (elem)aux;
        	}
        }
        // load b
        for (int i=0; i<N_HOR; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
        	b[i] = (elem)aux;
        }
        // load expectedResult
        for (int i=0; i<N_HOR; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            expectedResult[i] = aux;
        }
        cholHW(A, b, z);

        // Special care needed for zo (sometimes is all 0s and sometimes its the previous z)
        //cgradHW(A, b, zo, z, tol)
        //minresHW(A, b, zo, z, iterLS, tol);
        //cout << "sample number : " << sample << endl;
        for (int i=0; i<N_HOR; i++){
        	double error = fabs(expectedResult[i] -(double)z[i]);
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << z[i] << "\t";
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
