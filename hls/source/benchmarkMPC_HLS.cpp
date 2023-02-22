#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <iomanip>

#include "specs.h"
#include "tictoc.h"
#include "mmultSW.h"
//using namespace std;

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#include <CL/cl2.hpp>

#define DISPLAY
extern "C"{
void mpcHW(elem x_axi[N_SYS], elem r_axi[P_SYS], elem xmin_axi[N_SYS], elem xmax_axi[N_SYS], elem umin_axi[M_SYS], elem umax_axi[M_SYS], elem Acal_axi[N_SYS*N_HOR][N_SYS], elem AcalQOcal_axi[N_SYS][N_HOR], elem H_axi[N_HOR][N_HOR], elem M_hat_axi[I_AUX][N_HOR], elem L_invLast_axi[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u_axi[M_SYS]);
}

int main(int argc, char* argv[]) {

	std::cout << std::setprecision(15);
	std::cout << "MPC testbench" << std::endl;

	int iterPDIP = 20;
	int iterLS = 30;
	int errors = 0;
	double expectedU[M_SYS];
	double expectedX[N_SYS];

	elem A[N_SYS][N_SYS];
	elem B[N_SYS][M_SYS];
    elem x_next[N_SYS];
    elem x_nextAux[N_SYS];


	// threshold for the difference between the calculated results and the correct ones.
	double threshold = 1e-4;
	// tolerance for the linear solver in PDIP algorithm
	elem tol = 1e-9;
	// warmup samples to remove
	// warmup is removed for preliminary analysis, all the samples are
	// saved before MPCTimer.rmWarmup(warmupSamples);
	int warmupSamples = 100;

    //TARGET_DEVICE macro needs to be passed from gcc command line
    if(argc != 4) {
		std::cout << "Usage: " << argv[0] <<" <xclbin> <samples.bin> <times.csv>" << std::endl;
		return EXIT_FAILURE;
	}

    char* xclbinFilename = argv[1];


    std::vector<cl::Device> devices;
    cl::Device device;
    std::vector<cl::Platform> platforms;
    bool found_device = false;

    //traversing all Platforms To find Xilinx Platform and targeted
    //Device in Xilinx Platform
    cl::Platform::get(&platforms);
    for(size_t i = 0; (i < platforms.size() ) & (found_device == false) ;i++){
        cl::Platform platform = platforms[i];
        std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
        if ( platformName == "Xilinx"){
            devices.clear();
            platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
	    if (devices.size()){
		    device = devices[0];
		    found_device = true;
		    break;
	    }
        }
    }
    if (found_device == false){
       std::cout << "Error: Unable to find Target Device "
           << device.getInfo<CL_DEVICE_NAME>() << std::endl;
       return EXIT_FAILURE;
    }

    // Creating Context and Command Queue for selected device
    cl::Context context(device);
    cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE);

    // Load xclbin
    std::cout << "Loading: '" << xclbinFilename << "'\n";
    std::ifstream bin_file(xclbinFilename, std::ifstream::binary);
    bin_file.seekg (0, bin_file.end);
    unsigned nb = bin_file.tellg();
    bin_file.seekg (0, bin_file.beg);
    char *buf = new char [nb];
    bin_file.read(buf, nb);
    std::cout << "Done loading xclbin file\n";

    // Creating Program from Binary File
    std::cout << "Creating Program from Binary File\n";
    cl::Program::Binaries bins;
    bins.push_back({buf,nb});
    devices.resize(1);
    cl::Program program(context, devices, bins);
    std::cout << "DONE Creating Program from Binary File\n";

    // This call will get the kernel object from program. A kernel is an
    // OpenCL function that is executed on the FPGA.
    cl::Kernel krnl_mpc(program,"mpcHW");


// elem x[N_SYS], elem r[P_SYS], elem xmin[N_SYS], elem xmax[N_SYS], elem umin[M_SYS], elem umax[M_SYS], elem Acal[N_SYS*N_QP][N_SYS],
// elem AcalQOcal[N_SYS][N_QP], elem H[N_QP][N_QP], elem M_hat[N_SYS*N_QP*2+2*N_QP][N_QP], elem L_invLast[N_SYS+P_SYS], int iterPDIP, int iterLS, elem u[N_SYS])

    // These commands will allocate memory on the Device. The cl::Buffer objects can
    // be used to reference the memory locations on the device.
    std::cout << "Creating Buffers\n";
    cl::Buffer buffer_x(context, CL_MEM_READ_ONLY, N_SYS * sizeof(elem));
    cl::Buffer buffer_r(context, CL_MEM_READ_ONLY, P_SYS * sizeof(elem));

    cl::Buffer buffer_xmin(context, CL_MEM_READ_ONLY, N_SYS * sizeof(elem));
    cl::Buffer buffer_xmax(context, CL_MEM_READ_ONLY, N_SYS * sizeof(elem));
    cl::Buffer buffer_umin(context, CL_MEM_READ_ONLY, M_SYS * sizeof(elem));
    cl::Buffer buffer_umax(context, CL_MEM_READ_ONLY, M_SYS * sizeof(elem));

    cl::Buffer buffer_Acal(context, CL_MEM_READ_ONLY, N_SYS*N_HOR*N_SYS * sizeof(elem));
    cl::Buffer buffer_AcalQOcal(context, CL_MEM_READ_ONLY, N_SYS*N_HOR * sizeof(elem));
    cl::Buffer buffer_L_invLast(context, CL_MEM_READ_ONLY, (N_SYS+P_SYS) * sizeof(elem));

    cl::Buffer buffer_H(context, CL_MEM_READ_ONLY, N_HOR*N_HOR * sizeof(elem));
    cl::Buffer buffer_M_hat(context, CL_MEM_READ_ONLY, I_AUX*N_HOR* sizeof(elem));

    //cl::Buffer buffer_iterPDIP(context, CL_MEM_READ_ONLY, sizeof(int));
    //cl::Buffer buffer_iterLS(context, CL_MEM_READ_ONLY, sizeof(int));

    cl::Buffer buffer_u(context, CL_MEM_WRITE_ONLY, M_SYS * sizeof(elem));

    //set the kernel Arguments
    std::cout << "Setting kernel arguments\n";
    int narg=0;
    krnl_mpc.setArg(narg++,buffer_x);
    krnl_mpc.setArg(narg++,buffer_r);

    krnl_mpc.setArg(narg++,buffer_xmin);
    krnl_mpc.setArg(narg++,buffer_xmax);
    krnl_mpc.setArg(narg++,buffer_umin);
    krnl_mpc.setArg(narg++,buffer_umax);

    krnl_mpc.setArg(narg++,buffer_Acal);
    krnl_mpc.setArg(narg++,buffer_AcalQOcal);


    krnl_mpc.setArg(narg++,buffer_H);
    krnl_mpc.setArg(narg++,buffer_M_hat);

    krnl_mpc.setArg(narg++,buffer_L_invLast);

    krnl_mpc.setArg(narg++,iterPDIP);
    krnl_mpc.setArg(narg++,iterLS);

    krnl_mpc.setArg(narg++,buffer_u);

    std::cout << "Mapping Buffers\n";
    //We then need to map our OpenCL buffers to get the pointers
    std::cout << "Buffer: x\n";
	elem *x= (elem *) q.enqueueMapBuffer (buffer_x , CL_TRUE , CL_MAP_WRITE , 0, N_SYS * sizeof(elem));

	std::cout << "Buffer: r\n";
	elem *r = (elem *) q.enqueueMapBuffer (buffer_r , CL_TRUE , CL_MAP_WRITE , 0, P_SYS * sizeof(elem));

	std::cout << "Buffer: xmin\n";
	elem *xmin = (elem *) q.enqueueMapBuffer (buffer_xmin , CL_TRUE , CL_MAP_WRITE , 0, N_SYS * sizeof(elem));

	std::cout << "Buffer: xmax\n";
	elem *xmax = (elem *) q.enqueueMapBuffer (buffer_xmax , CL_TRUE , CL_MAP_WRITE , 0, N_SYS * sizeof(elem));

	std::cout << "Buffer: umin\n";
	elem *umin = (elem *) q.enqueueMapBuffer (buffer_umin , CL_TRUE , CL_MAP_WRITE , 0, M_SYS * sizeof(elem));

	std::cout << "Buffer: umax\n";
	elem *umax = (elem *) q.enqueueMapBuffer (buffer_umax , CL_TRUE , CL_MAP_WRITE , 0, M_SYS * sizeof(elem));

	std::cout << "Buffer: Acal\n";
	elem *Acal = (elem *) q.enqueueMapBuffer (buffer_Acal , CL_TRUE , CL_MAP_WRITE , 0, N_SYS*N_HOR*N_SYS * sizeof(elem));

	std::cout << "Buffer: AcalQCal\n";
	elem *AcalQOcal = (elem *) q.enqueueMapBuffer (buffer_AcalQOcal , CL_TRUE , CL_MAP_WRITE , 0, N_SYS*N_HOR * sizeof(elem));

	std::cout << "Buffer: buffer_L_invLast\n";
	elem *L_invLast= (elem *) q.enqueueMapBuffer (buffer_L_invLast , CL_TRUE , CL_MAP_WRITE , 0, (N_SYS+P_SYS) * sizeof(elem));

	std::cout << "Buffer: H\n";
	elem *H = (elem *) q.enqueueMapBuffer (buffer_H , CL_TRUE , CL_MAP_WRITE , 0, N_HOR*N_HOR * sizeof(elem));

	std::cout << "Buffer: M_hat\n";
	elem *M_hat = (elem *) q.enqueueMapBuffer (buffer_M_hat , CL_TRUE , CL_MAP_WRITE , 0, I_AUX*N_HOR * sizeof(elem));



	std::cout << "Buffer: u\n";
	elem *u = (elem *) q.enqueueMapBuffer (buffer_u , CL_TRUE , CL_TRUE , 0,  M_SYS*sizeof(elem));

	std::cout << "Buffers mapped\n";



	std::ifstream samples(argv[2], std::ios::binary);
    //check to see that the file was opened correctly:
    if (!samples.is_open()) {
    	std::cerr << "There was a problem opening the input file: ";
    	std::cerr << argv[2] << std::endl;
        return 1;//exit or do additional error checking
    }

    std::ofstream csvTimes(argv[3], std::ios::app);
    //check to see that the file was opened correctly:
    if (!csvTimes.is_open()) {
        std::cerr << "There was a problem opening the output file: ";
        std::cerr << argv[2] <<  std::endl;
        return 1;//exit or do additional error checking
    }



    int nSamples = 0;
    int samplesI = 0;
    int samplesN = 0;


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
	tol = (elem)aux;

    // check if N and I match with N_QP and I_QP from specs.h
    if (N_HOR != samplesN){
    	std::cerr << "N_HOR from " << argv[1] << " (" << samplesN << ") does not match N_QP from specs.h (";
    	std::cerr << N_HOR << ")" << std::endl;
        return 1;
    }
    if (I_AUX != samplesI){
    	std::cerr << "I_AUX from " << argv[1] << " (" << samplesI << ")does not match I_QP from specs.h (";
    	std::cerr << I_AUX << ")" << std::endl;
        return 1;
    }

    // load xmin
    for (int r=0; r<N_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	xmin[r] = (elem)aux;
    	//std::cout << "xmin["<< r << "]: " << xmin[r] << std::endl;
    }
    // load xmax
    for (int r=0; r<N_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	xmax[r] = (elem)aux;
    	//std::cout << "xmax["<< r << "]: " << xmax[r] << std::endl;
    }
    // load umin
    for (int r=0; r<M_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	umin[r] = (elem)aux;
    	//std::cout << "umin["<< r << "]: " << umin[r] << std::endl;
    }
    // load umax
    for (int r=0; r<M_SYS; r++){
    	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    	umax[r] = (elem)aux;
    	//std::cout << "umax["<< r << "]: " << umax[r] << std::endl;
    }
    // load Acal
    for (int r=0; r<N_SYS*N_HOR; r++){
    	for (int c=0; c<N_SYS; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		Acal[r*N_SYS + c] = (elem)aux;
    		//std::cout << "Acal["<< r*N_SYS + c << "]: " << Acal[r*N_SYS + c] << std::endl;
    	}
    }
    // load AcalOmgOcal
    for (int r=0; r<N_SYS; r++){
    	for (int c=0; c<N_HOR; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		AcalQOcal[r*N_HOR + c] = (elem)aux;
    		//std::cout << "AcalQOcal["<< r*N_HOR + c << "]: " << AcalQOcal[r*N_HOR + c] << std::endl;
    	}
    }
    // load H
    for (int r=0; r<N_HOR; r++){
    	for (int c=0; c<N_HOR; c++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		H[r*N_HOR + c] = (elem)aux;
    		//std::cout << "H["<< r*N_HOR + c << "]: " << H[r*N_HOR + c] << std::endl;
    	}
    }
    // load M_hat
	for (int r=0; r<I_AUX; r++){
		for (int c=0; c<N_HOR; c++){
			samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
			M_hat[r*N_HOR + c] = (elem)aux;
			//std::cout << "M_hat["<< r*N_HOR + c << "]: " << M_hat[r*N_HOR + c] << std::endl;
		}
    }
    // load L_invLast
	for (int r=0; r<N_SYS+P_SYS; r++){
		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
		L_invLast[r] = (elem)aux;
		//std::cout << "L_invLast["<< r << "]: " << L_invLast[r] << std::endl;
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
	// load initial state
	for (int i=0; i<N_SYS; i++){
		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
		x[i] = (elem)aux;
	}

    // Data will be migrated to kernel space
    q.enqueueMigrateMemObjects({buffer_xmin,buffer_xmax,buffer_umin,buffer_umax},0/* 0 means from host*/);
    q.enqueueMigrateMemObjects({buffer_Acal,buffer_AcalQOcal,buffer_L_invLast},0/* 0 means from host*/);
    q.enqueueMigrateMemObjects({buffer_H,buffer_M_hat},0/* 0 means from host*/);

    //Launch the Kernel
    q.enqueueTask(krnl_mpc);

    // The result of the previous kernel execution will need to be retrieved in
    // order to view the results. This call will transfer the data from FPGA to
    // source_results vector
    q.enqueueMigrateMemObjects({buffer_u},CL_MIGRATE_MEM_OBJECT_HOST);

    q.finish();

    // tic toc timer
    tictoc MPCTimer(nSamples);

    for (int sample=0; sample<nSamples; sample++){

        // load r
        for (int i=0; i<P_SYS; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            r[i] = (elem)aux;
        }
        // load expectedU
        for (int i=0; i<P_SYS; i++){
        	samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
            expectedU[i] = aux;
        }
        // load expected x
    	for (int i=0; i<N_SYS; i++){
    		samples.read(reinterpret_cast<char*>(&aux), sizeof(double));
    		expectedX[i] = (elem)aux;
    	}
        MPCTimer.tic();
        // Load x and r
        q.enqueueMigrateMemObjects({buffer_x,buffer_r},0/* 0 means from host*/);
        // Launch the Kernel
        q.enqueueTask(krnl_mpc);
        // Read u
        q.enqueueMigrateMemObjects({buffer_u},CL_MIGRATE_MEM_OBJECT_HOST);
        // Wait for the kernel to finish
        q.finish();
        MPCTimer.toc();

        // X = A*X + B*U
        mmultSW (&A[0][0],x, x_next, N_SYS, N_SYS, 1);
        mmultSW (&B[0][0],u, x_nextAux, N_SYS, M_SYS, 1);
        for (int i=0; i<N_SYS; i++){
			x_next[i] += x_nextAux[i];
		}

        // exit if end of file
        if (samples.peek() == EOF) break;

        for (int i=0; i<N_SYS; i++){
            double error = fabs(expectedX[i] -(double)x_next[i]);
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                std::cout << "sample number : " << sample << " " << i << std::endl;
                std::cout << expectedX[i] << "\t";
                std::cout << x_next[i] << "\t";
                std::cout << error;
                std::cout << "\n";
#endif

            }
        }
        for (int i=0; i<N_SYS; i++){
            x[i] = x_next[i];
        }

    }

    samples.close();

    MPCTimer.exportTocs(csvTimes);
    csvTimes.close();

    // print important data
    std::cout << "Finished processing " << nSamples << " samples" << std::endl;
    std::cout << "Number of differences between expected and calculated " << errors << std::endl;
    std::cout << "Threshold: " << threshold << std::endl;
    std::cout << std::endl;

    MPCTimer.rmWarmup(warmupSamples);

    std::cout << "Mean time running application in software: \t\t" << MPCTimer.mean() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Standard deviation running application in software: \t" << MPCTimer.stdev() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Max running application in software: \t\t\t" << MPCTimer.max() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Min running application in software: \t\t\t" << MPCTimer.min() / 1e3 << "\xC2\xB5s" << std::endl;



    q.enqueueUnmapMemObject(buffer_xmin , xmin);
    q.enqueueUnmapMemObject(buffer_xmax , xmax);
    q.enqueueUnmapMemObject(buffer_umin , umin);
    q.enqueueUnmapMemObject(buffer_umax , umax);
    q.enqueueUnmapMemObject(buffer_Acal , Acal);
    q.enqueueUnmapMemObject(buffer_AcalQOcal , AcalQOcal);
    q.enqueueUnmapMemObject(buffer_L_invLast , L_invLast);
    q.enqueueUnmapMemObject(buffer_H , H);
    q.enqueueUnmapMemObject(buffer_M_hat , M_hat);
    q.enqueueUnmapMemObject(buffer_Acal , Acal);
    q.enqueueUnmapMemObject(buffer_x , x);
    q.enqueueUnmapMemObject(buffer_r , r);
    q.enqueueUnmapMemObject(buffer_u , u);
    q.finish();


    return errors;

}
