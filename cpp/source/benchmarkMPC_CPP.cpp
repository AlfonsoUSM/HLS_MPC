#include <iostream>
#include <fstream>
#include "specs.h"
#include <math.h>
#include "mpcSW.h"
#include "tictoc.h"
#include <iomanip>




using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "MPC testbench" << endl;

	elem x[N_SYS];
    elem x_next[N_SYS];
    elem x_nextAux[N_SYS];
	elem r[P_SYS];
	elem xmin[N_SYS];
	elem xmax[N_SYS];
	elem umin[P_SYS];
	elem umax[P_SYS];
	elem Acal[N_SYS*N_HOR][N_SYS];
	elem AcalQOcal[N_SYS][N_HOR];
	elem H[N_HOR][N_HOR];
	elem M_hat[I_AUX][N_HOR];
	elem L_invLast[(N_SYS+P_SYS)];
	elem u[M_SYS];
    elem A[N_SYS][N_SYS];
    elem B[N_SYS][M_SYS];
    double expectedResult[P_SYS];
	int iterPDIP = 0;
	int iterLS = 0;
	int errors = 0;
    // threshold for the difference between the calculated results and the correct ones.
    double threshold = 1e-4;
    // tolerance for the linear solver in PDIP algorithm
    elem tol = 1e-9;
    // warmup samples to remove
    int warmupSamples = 100;

    if (argc!=3){
        cerr << "Must specify .bin\n";
        cerr << "and/or\n";
        cerr << "Must specify .csv\n";
        return 1;
    }

    ifstream samples(argv[1], ios::binary);
    //check to see that the file was opened correctly:
    if (!samples.is_open()) {
        cerr << "There was a problem opening the input file: ";
        cerr << argv[1] << endl;
        return 1;//exit or do additional error checking
    }

    ofstream csvTimes(argv[2], ios::app);
    //check to see that the file was opened correctly:
    if (!csvTimes.is_open()) {
        cerr << "There was a problem opening the output file: ";
        cerr << argv[2] << endl;
        return 1;//exit or do additional error checking
    }



    double samplesN = 0;
    double samplesI = 0;
    double nSamples = 0;



    
    // get number N, I and number of samples
    samples.read(reinterpret_cast<char*>(&samplesN), sizeof(double));
    samples.read(reinterpret_cast<char*>(&samplesI), sizeof(double));
    samples.read(reinterpret_cast<char*>(&nSamples), sizeof(double));

    // check if N and I match with N_HOR and I_AUX from specs.h
    if (N_HOR != (int)samplesN){
        cerr << "N_HOR from " << argv[1] << " (" << samplesN << ") does not match N_HOR from specs.h (";
        cerr << N_HOR << ")" << endl;
        return 1;
    }
    if (I_AUX != (int)samplesI){
        cerr << "I_AUX from " << argv[1] << " (" << samplesI << ")does not match I_AUX from specs.h (";
        cerr << I_AUX << ")" << endl;
        return 1;
    }

    double aux;
    //load iterPDIP
    samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    iterPDIP = (int)aux;

    //load iterLS
    samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    iterLS = (int)aux;

    //load tol
    samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    tol = (elem)aux;

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
    for (int r=0; r<P_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
        umin[r] = (elem)aux;
    }
    // load umax
    for (int r=0; r<P_SYS; r++){
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
    // load A
	for (int r=0; r<N_SYS; r++){
		for (int c=0; c<N_SYS; c++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			A[r][c] = (elem)aux;
		}
    }
    // load B
	for (int r=0; r<N_SYS; r++){
		for (int c=0; c<M_SYS; c++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			B[r][c] = (elem)aux;
		}
    }

    // tic toc timer
    tictoc MPCTimer(nSamples);

    for (int sample=0; sample<nSamples; sample++){
    	// exit if end of file
        
        if (samples.peek() == EOF) break;
        // load x
        for (int i=0; i<N_SYS; i++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            x[i] = (elem)aux;
        }
        
        // load r file
        for (int i=0; i<P_SYS; i++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            r[i] = (elem)aux;
        }
        // load expectedResult
        for (int i=0; i<P_SYS; i++){
            samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            expectedResult[i] = (elem)aux;
        }

        if (sample>0){
            
            for (int i=0; i<N_SYS; i++){
                double error = fabs(x[i] -(double)x_next[i]);
                if ((error > threshold) | (error != error)){
                    errors++;
#ifdef DISPLAY
                    cout << "sample number : " << sample << " " << i << endl;
                    cout << x[i] << "\t";
                    cout << x_next[i] << "\t";
                    cout << error;
                    cout << "\n";
#endif

                }
            }
            for (int i=0; i<N_SYS; i++){
                x[i] = x_next[i];
            }
        }
        

        MPCTimer.tic();
        mpcSW(x, r, xmin, xmax, umin, umax, Acal, AcalQOcal, H, M_hat, L_invLast, iterPDIP, iterLS, tol, u);
        MPCTimer.toc();
        
        // X = A*X + B*U
        mmultSW (&A[0][0],x, x_next, N_SYS, N_SYS, 1);
        mmultSW (&B[0][0],u, x_nextAux, N_SYS, M_SYS, 1);
        for (int i=0; i<N_SYS; i++){
            x_next[i] += x_nextAux[i];
        }
        
        //cout << "sample number : " << sample << endl;


    }

    samples.close();

    // print important data
    cout << "Finished processing " << nSamples << " samples" << endl;
    cout << "Number of differences between expected and calculated:\t" << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;
    
    // save times with warmup
    MPCTimer.exportTocs(csvTimes);
    MPCTimer.rmWarmup(warmupSamples);

    std::cout << "Mean time processing a sample:\t\t" << MPCTimer.mean() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Standard deviation processing a sample:\t" << MPCTimer.stdev() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Max running processing a sample:\t" << MPCTimer.max() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Min running processing a sample:\t" << MPCTimer.min() / 1e3 << "\xC2\xB5s" << std::endl;
    
    
    csvTimes.close();
    return errors;
}
