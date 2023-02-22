#ifndef TICTOC_H
#define TICTOC_H
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#define PC

#include "specs.h"
#ifdef PC
#include <chrono>
using namespace std::chrono;
typedef high_resolution_clock::time_point tic_dt;
#else
#include "sds_lib.h"
typedef uint64_t tic_dt;
#endif



#ifdef PC

    class tictoc
    {

    private:

    public:
        std::vector<double> tocs;
        tic_dt start;
        tic_dt end;
        tictoc(int len);
        void reset();
        inline void tic(){
            start = high_resolution_clock::now();
        }
        inline void toc(){
            end = high_resolution_clock::now();
            double interval = duration_cast<nanoseconds>(end - start).count();
            tocs.push_back(interval);
        }
        void exportTocs(std::ofstream &csvOut);
        void rmWarmup(int n);
        double mean();
        double stdev();
        double max();
        double min();
    };
#else

    class tictoc {

    private:

    public:
        std::vector<double> tocs;
        tic_dt start;
        tic_dt end;
        double freq = sds_clock_frequency();
        int pos;
        tictoc(int len);
        void reset();
        inline void tic(){
            start = sds_clock_counter();
        }
        inline void toc(){
            end = sds_clock_counter();
            tocs.push_back((end - start)/freq);
            pos++;
        }
        uint64_t mean();
        double stdev();
        uint64_t max();
        uint64_t min();
    };
#endif
#endif      // TICTOC_H
