/**
 * @file    main.cu
 * @author  Yohai Meiron <ymeiron@pku.edu.cn>
 * @version 1.0
 */

// If you have CUDA then compile with:
// nvcc main.cu -lm -O3 -arch=sm_20 -o output
// Otherwise enable OpenMP and compile with GCC:
// g++ -x c++ -O3 -o output main.cu -DOMP -fopenmp -lgomp -I/home/ym
// The -I is the path to the parent directory where thrust is.

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>



#define CUDA

#ifdef OMP
    #error Sorry, OpenMP is currently disabled.
    #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_BACKEND_OMP
    #undef CUDA
    #define PARALLEL_GET_TID omp_get_thread_num()
    #define PARALLEL_ADVANCE omp_get_num_threads()
    #define __global__
//     #include <complex>
#else
    #define PARALLEL_GET_TID threadIdx.x + blockIdx.x * blockDim.x
    #define PARALLEL_ADVANCE blockDim.x * gridDim.x
//     #include "cuda_complex.hpp"
#endif
// #define complex complex<Real>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
//#include <thrust/remove.h>
#include <thrust/partition.h>

#include "common.hpp"
#include "io.hpp"
#include "integrate.hpp"

#ifdef MEX
    #define method mex
    #include "mex.hpp"
#elif defined(SCF)
    #define method scf
    #include "scf.hpp"
#endif

using namespace std;



/////////////////////////////////////// DBG
#include <sys/timeb.h>
int AccurateGetTime() {
    timeb tb;
    ftime(&tb);
    return tb.millitm + (tb.time & 0xfffff) * 1000;
}

int AccurateTimeDiff(int End, int Start) {
    int Span = End - Start;
    if (Span < 0) Span += 0x100000 * 1000;
    return Span;
}
int DBG_EVENT_START;
int DBG_EVENT_END;
int DBG_EVENT_CALLS = 0;


//#define DBG_SUPPRESS_OUTPUT

#define DBG_MACRO_EXIT {exit(1);}
#define DBG_MACRO_SKIP {return 0;}
#ifndef DBG_SUPPRESS_OUTPUT
#define DBG_CERR cerr
#define DBG_STDERR stderr
#define DBG_MACRO_INIT {                                                       \
    DBG_EVENT_CALLS++;                                                         \
    fprintf(DBG_STDERR, "(%d) Initializing CPU clock\n", DBG_EVENT_CALLS);     \
    DBG_EVENT_START = AccurateGetTime();                                       \
}
#define DBG_MACRO_MSG(MSG) {                                                   \
    DBG_EVENT_END = AccurateGetTime();                                         \
    fprintf(DBG_STDERR, "(%d) Previous event time: %d ms\n", DBG_EVENT_CALLS,  \
      AccurateTimeDiff(DBG_EVENT_END, DBG_EVENT_START));                       \
    DBG_CERR << MSG << endl;                                                   \
}
cudaEvent_t DBG_CUDA_EVENT_START, DBG_CUDA_EVENT_END;
float DBG_CUDA_ELAPSED_TIME;
#define DBG_MACRO_CUDA_INIT {                                                  \
    DBG_EVENT_CALLS++;                                                         \
    fprintf(DBG_STDERR, "(%d) Initializing GPU clock\n", DBG_EVENT_CALLS);     \
    cudaEventCreate(&DBG_CUDA_EVENT_START);                                    \
    cudaEventCreate(&DBG_CUDA_EVENT_END);                                      \
    cudaEventRecord(DBG_CUDA_EVENT_START, 0);                                  \
}
#define DBG_MACRO_CUDA_MSG(MSG) {                                              \
    cudaEventRecord(DBG_CUDA_EVENT_END, 0);                                    \
    cudaEventSynchronize(DBG_CUDA_EVENT_END);                                  \
    cudaEventElapsedTime(&DBG_CUDA_ELAPSED_TIME, DBG_CUDA_EVENT_START,         \
      DBG_CUDA_EVENT_END);                                                     \
    fprintf(DBG_STDERR, "(%d) Previous event time: %.2f ms\n", DBG_EVENT_CALLS,\
      DBG_CUDA_ELAPSED_TIME);                                                  \
    DBG_CERR << MSG << endl;                                                   \
}
#else
ostream DBG_CERR(0);
FILE *DBG_STDERR = fopen("/dev/null", "w");
#define DBG_MACRO_INIT
#define DBG_MACRO_MSG(MSG)
#define DBG_MACRO_CUDA_INIT
#define DBG_MACRO_CUDA_MSG(MSG)
#endif
/////////////////////////////////////// DBG

// GLOBAL VARIABLES
extern Real ConstantStep;
extern Real T, Step, dT1, dT2, Tcrit;
extern int N;

extern Particle *hostP;
extern thrust::device_vector<Particle> P;

extern thrust::device_vector<vec3>   F0;
extern thrust::device_vector<Real>   Potential;
extern thrust::device_vector<vec3>   F1;

 
int NSteps = 0, SnapNumber = 0;
//TODO Not use global variables... They cause an 'unload of CUDA runtime failed'
// error message (which does not actually affect performance).


struct ReorderingFunctor {
    __host__ __device__ bool operator() (const Particle &lhs, const Particle &rhs) {
        return (lhs.ID <= rhs.ID);
    }
};

Real CalculateStepSize() {
    return ConstantStep;
}

struct KineticEnergyFunctor {
    __host__ __device__ Real operator() (const Particle &p) const {return 0.5*p.m*p.vel.abs2();}
};

Real KineticEnergy() {
    return thrust::transform_reduce(
      P.begin(),
      P.end(),
      KineticEnergyFunctor(),
      (Real)0, // It must be clear to the function that this zero is a Real.
      thrust::plus<Real>()
    );
}

void DisplayInformation(thrust::device_vector<Particle> *P, Real *Potential) {
    Real Ek = KineticEnergy();
    Real Ep = method::PotentialEnergy();
    Real Energy = Ek + Ep;
    printf(" TIME =%6.2f  NSTEPS =%6d  ENERGY =%20.16f   N = %d\n", T, NSteps, Energy, (*P).size());
    fflush(stdout);

}

void PrepareSnapshot(thrust::device_vector<Particle> *P, Particle *hostP) {
    Particle *P_ptr = thrust::raw_pointer_cast((*P).data());
    cudaMemcpy(hostP, P_ptr, N * sizeof(Particle), cudaMemcpyDeviceToHost);
    thrust::sort(hostP, hostP+N, ReorderingFunctor());
}

#define PTR(x) (thrust::raw_pointer_cast((x).data()))
#define numberoftries 10

int main(int argc, char *argv[]) {
    cerr << "Welcome to ETICS..." << endl;
#ifdef MEX
    cerr << "Using method: MEX" << endl;
    cerr << "LMAX=" << LMAX << endl;
#elif defined(SCF)
    cerr << "Using method: SCF" << endl;
    cerr << "LMAX=" << LMAX << endl;
    cerr << "NMAX=" << NMAX << endl;
#endif

    string FileName;
    int DeviceID = 0;

    ParametersStruct Params;
    ParseInput(argc, argv, &Params);
    N = Params.N;
    FileName =Params.FileName;
    Tcrit =Params.Tcrit;
    ConstantStep =Params.ConstantStep;
    DeviceID =Params.DeviceID;
    dT1 = Params.dT1;
    dT2 = Params.dT2;

    if (DeviceID >= 0) {
        if (cudaSetDevice(DeviceID) != cudaSuccess) {
            cerr <<  "Problem opening device (ID=" << DeviceID << ")" << endl;
            exit(1);
        }
    } else {
        cerr << "Skipping call to cudaSetDevice." << endl;
    }

     method::Init(N, 180, 64, 2605, 384);

    // Read an input file and initialize the global particle structure.
    ReadICs(FileName, N, Params.Skip, &hostP);
    InitilizeIntegratorMemory();
    CommitParticles();


//     for (int i=0; i<10; i++) method::CalculateGravity(PTR(P), N, PTR(Potential), PTR(F0)); // heat up the GPU
// //     K3Search(N);
// //     exit(4);
// 
//     double total=0;
// 
//     for (int i=0; i<numberoftries; i++) {
//         float milliseconds = 0;
//         cudaEvent_t start, stop;
//         cudaEventCreate(&start);
//         cudaEventCreate(&stop);
//         cudaEventRecord(start);
//         method::CalculateGravity(PTR(P), N, PTR(Potential), PTR(F0));
//         cudaEventRecord(stop);
//         cudaEventSynchronize(stop);
//         cudaEventElapsedTime(&milliseconds, start, stop);
//         total += milliseconds;
//     }
//     total /= numberoftries;
// 
//         thrust::host_vector<vec3>   F0_host(F0);
//         thrust::host_vector<double> PotHost(Potential);
//     for (int i=0; i<N; i+=1000) {
//         printf("%06d %+.16e %+.16e %+.16e %+.16e\n", i, F0_host[i].x, F0_host[i].y, F0_host[i].z, PotHost[i]);
//     }
//     printf("Time: %f\n", total);
//     
// //     OptimizeKernel4LaunchConfig();
//     
//     exit(0);
//     


    // Calculate the forces from the initial data.
    method::CalculateGravity(PTR(P), N, PTR(Potential), PTR(F0));
    CommitForces();

    // More initializations.
    Real NextOutput = 0, NextSnapshot = 0;
    T  = 0;
    DBG_MACRO_INIT
    Step = CalculateStepSize();
    DBG_MACRO_MSG("Finished [first] call to CalculateStepSize().")

    while (T <= Tcrit) {
        if (T >= NextOutput) {
            DisplayInformation(&P, PTR(Potential));
            NextOutput += dT1;
        }
        if (T >= NextSnapshot) {
#ifdef MEX
            PrepareSnapshot(&P, hostP);
#endif
            WriteSnapshot(Params.Prefix, SnapNumber, hostP, N, T);
            SnapNumber++;
            NextSnapshot += dT2;
        }

        // Take the drift step.
        DriftStep();

        // Calculate the forces in the new positions.
        method::CalculateGravity(PTR(P), N, PTR(Potential), PTR(F1));

        // Finish by taking the kick step.
        // The kick functor also "commits" the predicted forces into the "acc" member.
        KickStep();

        // N particles were implicitly propagated in this iteration.
        NSteps += 1;

        // Advance global time.
        T += Step;

        // Calculate the next step.
        Step = CalculateStepSize();
    }
    return 0;
}
