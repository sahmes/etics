#pragma once
#include "common.hpp"
#include "mathaux.hpp"
namespace etics {
    namespace scf {
        __global__ void LoadParticlesToCache(Particle *P, int N, vec3 ExpansionCenter=vec3(0,0,0));
        __global__ void CalculatePhi0l(int l);
        __global__ void CalculateCoefficientsPartial(int n, int l, Complex *PartialSum);
        template<int Mode> __device__ void CalculateGravityTemplate(int i, Complex *A, vec3 *F, Real *Potential);
        __global__ void CalculateGravityFromCoefficients(Real *Potential, vec3 *F);
        void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
        void Init(int N, int k3gs_new=0, int k3bs_new=0, int k4gs_new=0, int k4bs_new=0); // to be removed!
        void GlobalInit();
        void GuessLaunchConfiguration(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new);

        struct CacheStruct {
            int N;
            Real *xi;
            Real *Phi0l;
            Real *Wprev1;
            Real *Wprev2;
            Real *costheta;
            Real *sintheta_I;
            Complex *Exponent;
            Real *mass;
        };

        class scfclass {
          public:
            scfclass();
            ~scfclass();
            void CalculateCoefficients();
            void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
            void Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
          private:
            Complex A_h[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
            CacheStruct Cache_h;
            Complex *PartialSum;
            Complex *PartialSum_h;
            int k3gs, k3bs, k4gs, k4bs;
            void GetGpuLock();
            void ReleaseGpuLock();
            void SendCachePointersToGPU();
            void SendCoeffsToGPU();
            void CalculateCoefficients(int n, int l);
        };
    }
}
