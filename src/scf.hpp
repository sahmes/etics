#pragma once
#include "common.hpp"
#include "mathaux.hpp"
namespace scf {
    void InitializeCache(int N);
    __global__ void LoadParticlesToCache(Particle *P, int N);
    __global__ void CalculatePhi0l(int l);
    __global__ void CalculateCoefficientsPartial(int n, int l, Complex *PartialSum);
    void CalculateCoefficients(int n, int l);
    void CalculateCoefficients();
    template<int Mode> __device__ void CalculateGravityTemplate(int i, Complex *A, vec3 *F, Real *Potential);
    __global__ void CalculateGravityFromCoefficients(Real *Potential, vec3 *F);
    void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
    void Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
    Real PotentialEnergy();
}
