#pragma once
#include "common.hpp"
#include "mathaux.hpp"
namespace etics {
    namespace scf {
        __global__ void LoadParticlesToCache(Particle *P, int N, vec3 Offset=vec3(0,0,0));
        __global__ void CalculatePhi0l(int N, int l);
        __global__ void CalculateCoefficientsPartial(int N, int n, int l, Complex *PartialSum);
        template<int Mode> __host__ __device__ void CalculateGravityTemplate(Real xi, Real costheta, Complex Exponent, Complex *A, vec3 *F, Real *Potential);
        __global__ void CalculateGravityFromCoefficients(int N, Real *Potential, vec3 *F);
        void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
        void Init(int Nmax, int k3gs_new=0, int k3bs_new=0, int k4gs_new=0, int k4bs_new=0); // to be removed!
        void GlobalInit();
        void GuessLaunchConfiguration(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new);

        struct CacheStruct {
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
            void SendCachePointersToGPU();
            void SendCoeffsToGPU();
            void LoadParticlesToCache(Particle *P, int N, vec3 Offset=vec3(0,0,0));
            void CalculateCoefficients();
            void CalculateGravityFromCoefficients(Real *Potential, vec3 *F);
            void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
            void CalculateGravityAtPoint(vec3 Pos, Real *Potential, vec3 *F);
            void Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
            void SetLaunchConfiguration(int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
            void GetCoefficients(Complex *A);
            void GetGpuLock();
            void ReleaseGpuLock();
//           private: // we'll make them private in the final version, now we need access for debugging
            int N;
            int Nmax;
            Complex A_h[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
            CacheStruct Cache_h;
            Complex *PartialSum;
            Complex *PartialSum_h;
            int k3gs, k3bs, k4gs, k4bs;
            void CalculateCoefficients(int n, int l);
        };
    }
}
