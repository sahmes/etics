#include "common.hpp"
#include "mathaux.hpp"
#include "scf.hpp"
#include "multicenter.hpp"
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>

etics::scf::scfclass SCFObject1;
etics::scf::scfclass SCFObject2;

vec3 *F_tmp;
vec3 *F_tmp2;
Real *Potential_tmp;
Real *Potential_tmp2;

void etics::multicenter::Init(int Nmax, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new) { // to be removed!!
    cudaMalloc((void**)&F_tmp, Nmax*sizeof(vec3));
    cudaMalloc((void**)&F_tmp2, Nmax*sizeof(vec3));
    cudaMalloc((void**)&Potential_tmp, Nmax*sizeof(Real));
    cudaMalloc((void**)&Potential_tmp2, Nmax*sizeof(Real));
    SCFObject1.Init(Nmax, k3gs_new, k3bs_new, k4gs_new, k4bs_new);
//     SCFObject1.OriginPos = vec3(5,0,0);
//     SCFObject1.OriginVel = vec3(0,0.1,0);

    SCFObject2.Init(Nmax, k3gs_new, k3bs_new, k4gs_new, k4bs_new);
}

void etics::multicenter::CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F) {
    SCFObject1.GetGpuLock();
    SCFObject1.SendCachePointersToGPU();
    SCFObject1.LoadParticlesToCache(P, N/2);
    SCFObject1.CalculateCoefficients();
    SCFObject1.ReleaseGpuLock();
    
    
    SCFObject2.GetGpuLock();
    SCFObject2.SendCachePointersToGPU();
    SCFObject2.LoadParticlesToCache(P+N/2, N/2);
    SCFObject2.CalculateCoefficients();
    SCFObject2.ReleaseGpuLock();
    
    
    SCFObject1.GetGpuLock();
    SCFObject1.SendCoeffsToGPU();
    SCFObject1.SendCachePointersToGPU();
    SCFObject1.LoadParticlesToCache(P, N);
    SCFObject1.CalculateGravityFromCoefficients(Potential_tmp, F_tmp);
    SCFObject1.ReleaseGpuLock();

    SCFObject2.GetGpuLock();
    SCFObject2.SendCoeffsToGPU();
    SCFObject2.SendCachePointersToGPU();
    SCFObject2.LoadParticlesToCache(P, N);
    SCFObject2.CalculateGravityFromCoefficients(Potential_tmp2, F_tmp2);
    SCFObject2.ReleaseGpuLock();
    
thrust::plus<Real> op;
thrust::transform(thrust::device, Potential_tmp, Potential_tmp+N, Potential_tmp2, Potential_tmp, thrust::plus<Real>());


    cudaMemcpy(Potential, Potential_tmp, N*sizeof(Real), cudaMemcpyDeviceToDevice);
    cudaMemcpy(F,         F_tmp,         N*sizeof(vec3), cudaMemcpyDeviceToDevice);
}