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
#include <mpi.h>


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
#include <thrust/inner_product.h>

#include "common.hpp"
#include "io.hpp"
#include "ic.hpp"
#include "integrate.hpp"

#ifdef MEX
    #define method mex
    #include "mex.hpp"
#elif defined(SCF)
    #define method scf
    #include "scf.hpp"
#endif

using namespace std;
using namespace etics;

// GLOBAL VARIABLES
int MyRank, NumProcs;
extern Real ConstantStep;
extern Real T, Step, dT1, dT2, Tcrit;
extern int N;

extern Particle *P_h;
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
      P.begin(), P.end(),
      KineticEnergyFunctor(),
      (Real)0, // It must be clear to the function that this zero is a Real.
      thrust::plus<Real>()
    );
}

struct PotentialEnergyFunctor {
    __host__ __device__ Real operator() (const Particle &p, const Real &Potential) const {return p.m*Potential;}
};

Real PotentialEnergy() {
    return 0.5*thrust::inner_product(
      P.begin(), P.end(),
      Potential.begin(),
      (Real)0,
      thrust::plus<Real>(),
      PotentialEnergyFunctor()
    );
}

void DisplayInformation(thrust::device_vector<Particle> *P, Real *Potential) {
    Real Ek = KineticEnergy();
    Real Ep = PotentialEnergy();
    Real Energy = Ek + Ep;

    Real TotalEnergy;
    MPI_Reduce(&Energy, &TotalEnergy, 1, MPI_ETICS_REAL, MPI_SUM, 0, MPI_COMM_WORLD);

    int N=(*P).size(), TotalN;
    MPI_Reduce(&N, &TotalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (MyRank==0) {
        printf(" TIME =%6.2f  NSTEPS =%6d  ENERGY =%20.16f   N = %d\n", T, NSteps, TotalEnergy, TotalN);
        fflush(stdout);
    }
}

#define PTR(x) (thrust::raw_pointer_cast((x).data()))

void PrepareSnapshot(thrust::device_vector<Particle> *P, Particle **ParticleList, int *CurrentTotalN) {
//     Particle *P_ptr = thrust::raw_pointer_cast((*P).data());
//     cudaMemcpy(ParticleList, P_ptr, N * sizeof(Particle), cudaMemcpyDeviceToHost);
    thrust::host_vector<Particle> LocalList((*P));
    int LocalBufferSize = LocalList.size()*sizeof(Particle);
    int BufferSizes[NumProcs];
    MPI_Gather(&LocalBufferSize, 1, MPI_INT, BufferSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int Displacements[NumProcs];
    int TotalN = 0;
    if (MyRank==0) {
        for (int p = 0; p < NumProcs; p++) TotalN += BufferSizes[p]/sizeof(Particle);
        Displacements[0] = 0;
        for (int p = 1; p < NumProcs; p++) Displacements[p] = Displacements[p-1] + BufferSizes[p-1];
        *ParticleList = new Particle[TotalN];
    }
    MPI_Gatherv(PTR(LocalList), LocalBufferSize, MPI_BYTE, *ParticleList, BufferSizes, Displacements, MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef MEX
    thrust::sort(*ParticleList, (*ParticleList)+TotalN, ReorderingFunctor());
#endif
    *CurrentTotalN = TotalN;
}


#define numberoftries 10

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
    MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);

    if (MyRank==0) {
        cerr << "Welcome to ETICS..." << endl;
#ifdef MEX
        cerr << "Using method: MEX" << endl;
        cerr << "LMAX=" << LMAX << endl;
#elif defined(SCF)
        cerr << "Using method: SCF" << endl;
        cerr << "LMAX=" << LMAX << endl;
        cerr << "NMAX=" << NMAX << endl;
#endif
    }

    string Filename;
    int DeviceID = 0;

    ParametersStruct Params;
    // Instead of reading the input file with MyRank=0 and broadcast the result, we let every rank read the file. This probably saves ~20 lines of ugly MPI code.
    ParseInput(argc, argv, &Params);
    N = Params.N; // total; will be divided by number of processes
    Filename = Params.Filename;
    Tcrit = Params.Tcrit;
    ConstantStep = Params.ConstantStep;
    DeviceID = Params.DeviceID;
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

    // Read an input file and initialize the global particle structure.
    Particle *FullList;
    if (MyRank==0) {
        if ((Filename.compare("_nofile_")==0) || (Filename.compare("_hernquist_")==0)) {
            cout << "Generating a Hernquist sphere..." << endl;
            etics::ic::hernquist(N, Params.Seed, &FullList);
            cout << "Done." << endl;
        } else if (Filename.compare("_plummer_")==0) {
            cout << "Generating a Plummer sphere..." << endl;
            etics::ic::plummer(N, Params.Seed, &FullList);
            cout << "Done." << endl;
        } else if (Filename.compare("_launch_config_")==0) {
            if (NumProcs > 1) {
                cerr << "Can only optimize launch configuration when a single MPI process is running!" << endl;
                exit(1);
            }
            int k3gs, k3bs, k4gs, k4bs;
            scf::OptimizeLaunchConfiguration(N, &k3gs, &k3bs, &k4gs, &k4bs);
            return 0;
        }
        else ReadICs(Filename, N, Params.Skip, &FullList);
    }

    int LocalN = N / NumProcs;
    int Remainder = N - LocalN*NumProcs;
    if (MyRank==NumProcs-1) LocalN += Remainder;
    P_h = new Particle[LocalN];
    int BufferSizes[NumProcs];
    int Displacements[NumProcs];
    if (MyRank==0) {
        for (int p = 0; p < NumProcs; p++) BufferSizes[p] = (N / NumProcs)*sizeof(Particle);
        BufferSizes[NumProcs-1] += Remainder*sizeof(Particle);
        Displacements[0] = 0;
        for (int p = 1; p < NumProcs; p++) Displacements[p] = Displacements[p-1] + BufferSizes[p-1];
    }
    MPI_Scatterv(FullList, BufferSizes, Displacements, MPI_BYTE, P_h, LocalN*sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD); 

    if (MyRank==0) free(FullList);
    N = LocalN;

    method::Init(N, 180, 64, 2605, 384);
#warning hardcoded launch configuration
    InitilizeIntegratorMemory();
    CommitParticles();

    // Calculate the forces from the initial data.
    method::CalculateGravity(PTR(P), N, PTR(Potential), PTR(F0));
    CommitForces();

    // More initializations.
    Real NextOutput = 0, NextSnapshot = 0;
    T  = 0;
    Step = CalculateStepSize();

    while (T <= Tcrit) {
        if (T >= NextOutput) {
            DisplayInformation(&P, PTR(Potential));
            NextOutput += dT1;
        }
        if (T >= NextSnapshot) {
            int CurrentTotalN;
            PrepareSnapshot(&P, &FullList, &CurrentTotalN);
            if (MyRank==0) {
                WriteSnapshot(Params.Prefix, SnapNumber, FullList, CurrentTotalN, T);
                free(FullList);
            }
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
    MPI_Finalize();
    return 0;
}
