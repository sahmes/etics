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
Real ConstantStep = 0.001953125;
Real T, Step, dT1, dT2, Tcrit, FileTime;
int NSteps = 0, FileSnapshotNum;

struct ReorderingFunctor {
    __host__ __device__ bool operator() (const Particle &lhs, const Particle &rhs) {
        return (lhs.ID <= rhs.ID);
    }
};

Real CalculateStepSize() {
    return ConstantStep;
}

void DisplayInformation(Integrator IntegratorObj) {
    Real Ek = IntegratorObj.KineticEnergy();
    Real Ep = IntegratorObj.PotentialEnergy();
    Real Energy = Ek + Ep;

    Real TotalEnergy;
    MPI_Reduce(&Energy, &TotalEnergy, 1, MPI_ETICS_REAL, MPI_SUM, 0, MPI_COMM_WORLD);

    int N=IntegratorObj.GetN(), TotalN;
    MPI_Reduce(&N, &TotalN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (MyRank==0) {
        printf(" TIME =%6.2f  NSTEPS =%6d  ENERGY =%20.16f   N = %d\n", T, NSteps, TotalEnergy, TotalN);
        fflush(stdout);
    }
}

void PrepareSnapshot(Integrator IntegratorObj, Particle **ParticleList, int *CurrentTotalN) {
    Particle *LocalList;
    int LocalBufferSize;
    IntegratorObj.CopyParticlesToHost(&LocalList, &LocalBufferSize);
    LocalBufferSize *= sizeof(Particle);
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
    MPI_Gatherv(LocalList, LocalBufferSize, MPI_BYTE, *ParticleList, BufferSizes,
 Displacements, MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef MEX
    thrust::sort(*ParticleList, (*ParticleList)+TotalN, ReorderingFunctor());
#endif
    *CurrentTotalN = TotalN;
}

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
    int N = Params.N; // total; will be divided by number of processes
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
        if ((Filename == "_nofile_") || (Filename == "_hernquist_")) {
            cout << "Generating a Hernquist sphere..." << endl;
            etics::ic::hernquist(N, Params.Seed, &FullList);
            FileSnapshotNum = 0;
            FileTime = 0;
            cout << "Done." << endl;
        } else if (Filename == "_plummer_") {
            cout << "Generating a Plummer sphere..." << endl;
            etics::ic::plummer(N, Params.Seed, &FullList);
            FileSnapshotNum = 0;
            FileTime = 0;
            cout << "Done." << endl;
        } else if (Filename == "_launch_config_") {
            if (NumProcs > 1) {
                cerr << "Can only optimize launch configuration when a single MPI process is running!" << endl;
                exit(1);
            }
            int k3gs, k3bs, k4gs, k4bs;
            scf::OptimizeLaunchConfiguration(N, &k3gs, &k3bs, &k4gs, &k4bs);
            return 0;
        }
        else {
            string InputFileSuffix = Filename.substr(Filename.find_last_of("."), Filename.length()-Filename.find_last_of("."));

            if ((InputFileSuffix==".h5part") || (InputFileSuffix==".hdf5") || (InputFileSuffix==".h5")) {
#ifndef ETICS_HDF5
                cerr << "Compiled without the \"ETICS_HDF5\" flag; cannot read input in this format." << endl;
                exit(1);
#else
                ReadICsHDF5(Filename, N, &FullList, &FileSnapshotNum, &FileTime);
#endif
            } else ReadICsASCII(Filename, N, &FullList, &FileSnapshotNum, &FileTime);
        }
    }
    
#ifndef ETICS_HDF5
    if (Params.Format == "hdf5") {
        cerr << "Compiled without the \"ETICS_HDF5\" flag; cannot output in requested format." << endl;
        exit(1);
    }
#endif
    if (!(Params.Format == "hdf5") && !(Params.Format == "ascii")) {
        cerr << "Requested output format unrecognized." << endl;
        exit(1);
    }

    int LocalN = N / NumProcs;
    int Remainder = N - LocalN*NumProcs;
    if (MyRank==NumProcs-1) LocalN += Remainder;
    Particle *LocalList = new Particle[LocalN];
    int BufferSizes[NumProcs];
    int Displacements[NumProcs];
    if (MyRank==0) {
        for (int p = 0; p < NumProcs; p++) BufferSizes[p] = (N / NumProcs)*sizeof(Particle);
        BufferSizes[NumProcs-1] += Remainder*sizeof(Particle);
        Displacements[0] = 0;
        for (int p = 1; p < NumProcs; p++) Displacements[p] = Displacements[p-1] + BufferSizes[p-1];
    }
    MPI_Scatterv(FullList, BufferSizes, Displacements, MPI_BYTE, LocalList, LocalN*sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (MyRank==0) free(FullList);
    N = LocalN;

    method::Init(N, 180, 64, 2605, 384);
#warning hardcoded launch configuration
    Integrator IntegratorObj(LocalList, N);

    // More initializations.
    Real NextOutput = 0, NextSnapshot = 0;
    T = FileTime;
    int SnapNumber = FileSnapshotNum;
    Step = CalculateStepSize();

    while (T <= Tcrit) {
        if (T >= NextOutput) {
            DisplayInformation(IntegratorObj);
            NextOutput += dT1;
        }
        if (T >= NextSnapshot) {
            int CurrentTotalN;
            PrepareSnapshot(IntegratorObj, &FullList, &CurrentTotalN);
            if (MyRank==0) {
                WriteSnapshotASCII(Params.Prefix, SnapNumber, FullList, CurrentTotalN, T);
                free(FullList);
            }
            SnapNumber++;
            NextSnapshot += dT2;
        }

        // Take the drift step.
        IntegratorObj.DriftStep(Step);

        // Calculate the forces in the new positions.
        IntegratorObj.CalculateGravity();

        // Finish by taking the kick step.
        // The kick functor also "commits" the predicted forces into the "acc" member.
        IntegratorObj.KickStep(Step);

        // N particles were implicitly propagated in this iteration.
        NSteps += 1;

        // Advance global time.
        T += Step;

        // Calculate the next step.
        Step = CalculateStepSize();
    }
    IntegratorObj.~Integrator();
    MPI_Finalize();
    return 0;
}
