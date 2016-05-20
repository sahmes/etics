#pragma once
#include "common.hpp"
#include "mathaux.hpp"
namespace etics {
    namespace scf {
        void InitializeCache(int N);
        void UpdateN(int N);
        __global__ void LoadParticlesToCache(Particle *P, int N, vec3 ExpansionCenter=vec3(0,0,0));
        __global__ void CalculatePhi0l(int l);
        __global__ void CalculateCoefficientsPartial(int n, int l, Complex *PartialSum);
        void CalculateCoefficients(int n, int l, Complex *A_h);
        void CalculateCoefficients(Complex *A_h);
        template<int Mode> __device__ void CalculateGravityTemplate(int i, Complex *A, vec3 *F, Real *Potential);
        __global__ void CalculateGravityFromCoefficients(Real *Potential, vec3 *F);
        void SendCoeffsToGPU(Complex *A_h);
        void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
        void Init(int N, int k3gs_new=0, int k3bs_new=0, int k4gs_new=0, int k4bs_new=0); // to be removed!
        void GlobalInit();
        void ZZZGuessLaunchConfigurationXXX(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new);

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
            void InitializeCache(int N);
            void CalculateCoefficients(int n, int l, Complex *A_h);
            void CalculateCoefficients(Complex *A_h);
            void SendCoeffsToGPU(Complex *A_h);
            void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
            void Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
            void GuessLaunchConfiguration(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new);
          private:
            Complex *PartialSum;
            Complex A_h[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
            CacheStruct Cache_h;
            Complex *PartialSum_h;
            int k3gs, k3bs, k4gs, k4bs;
            void GetGpuLock();
            void ReleaseGpuLock();
        };
    }
}
