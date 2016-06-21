#include "common.hpp"
#include "mathaux.hpp"
#include "scf.hpp"
#include "integrate.hpp"
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
Particle Origins[2];

void etics::multicenter::Init(int Nmax, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new) { // to be removed!!
    cudaMalloc((void**)&F_tmp, Nmax*sizeof(vec3));
    cudaMalloc((void**)&F_tmp2, Nmax*sizeof(vec3));
    cudaMalloc((void**)&Potential_tmp, Nmax*sizeof(Real));
    cudaMalloc((void**)&Potential_tmp2, Nmax*sizeof(Real));
    SCFObject1.Init(Nmax, k3gs_new, k3bs_new, k4gs_new, k4bs_new);
//     SCFObject1.OriginPos = vec3(5,0,0);
//     SCFObject1.OriginVel = vec3(0,0.1,0);

    SCFObject2.Init(Nmax, k3gs_new, k3bs_new, k4gs_new, k4bs_new);
    Particle p;
    p.m = 0;
    p.pos = vec3(-6.000000e-01 , 0 , 0);
    p.vel = vec3(0 , -3.803629e-03 , 0);
    Origins[0] = p;
    p.pos = vec3(2.940000e+01 , 0 , 0);
    p.vel = vec3(0 , 1.863778e-01 , 0);
    Origins[1] = p;
//     OriginsIntegrator = etics::Integrator(Origins, 2);
    
}

void etics::multicenter::CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F) {
    Particle PP[2];
    cudaMemcpy(PP, P+N-2, 2*sizeof(Particle), cudaMemcpyDeviceToHost);
//     std::cerr << PP[0].pos.x << std::endl;
//     std::cerr << PP[0].pos.y << std::endl;
//     std::cerr << PP[0].pos.z << std::endl;
//     std::cerr << PP[1].pos.x << std::endl;
//     std::cerr << PP[1].pos.y << std::endl;
//     std::cerr << PP[1].pos.z << std::endl;
    
    
//     #error all we did last time was to add the offset to the subroutines, now the energy is nan. must be some kind of self force
    
    SCFObject1.GetGpuLock();
    SCFObject1.SendCachePointersToGPU();
    SCFObject1.LoadParticlesToCache(P, 800000, PP[0].pos);
    SCFObject1.CalculateCoefficients();
    SCFObject1.ReleaseGpuLock();
    
    
    SCFObject2.GetGpuLock();
    SCFObject2.SendCachePointersToGPU();
    SCFObject2.LoadParticlesToCache(P+800000, 200000-2, PP[1].pos);
    SCFObject2.CalculateCoefficients();
    SCFObject2.ReleaseGpuLock();
    
    
    SCFObject1.GetGpuLock();
    SCFObject1.SendCoeffsToGPU();
    SCFObject1.SendCachePointersToGPU();
    SCFObject1.LoadParticlesToCache(P, N-2, PP[0].pos);
    SCFObject1.CalculateGravityFromCoefficients(Potential_tmp, F_tmp);
    SCFObject1.ReleaseGpuLock();
    
    SCFObject2.GetGpuLock();
    SCFObject2.SendCoeffsToGPU();
    SCFObject2.SendCachePointersToGPU();
    SCFObject2.LoadParticlesToCache(P, N-2, PP[1].pos);
    SCFObject2.CalculateGravityFromCoefficients(Potential_tmp2, F_tmp2);
    SCFObject2.ReleaseGpuLock();
    
    // now deal just with the central particles
    SCFObject1.GetGpuLock();
    SCFObject1.SendCoeffsToGPU();
    SCFObject1.SendCachePointersToGPU();
    SCFObject1.LoadParticlesToCache(P+N-1, 1, PP[0].pos); //instead of doing that, you need to "push" it at the end of the list
    SCFObject1.CalculateGravityFromCoefficients(Potential_tmp+N-1, F_tmp+N-1); // the last particle is the second center, so should feel the first's force
    SCFObject1.ReleaseGpuLock();

    SCFObject2.GetGpuLock();
    SCFObject2.SendCoeffsToGPU();
    SCFObject2.SendCachePointersToGPU();
    SCFObject2.LoadParticlesToCache(P+N-2, 1, PP[1].pos);
    SCFObject2.CalculateGravityFromCoefficients(Potential_tmp+N-2, F_tmp+N-2);
    SCFObject2.ReleaseGpuLock();
    
    
thrust::plus<Real> op;
thrust::transform(thrust::device, Potential_tmp, Potential_tmp+N-2, Potential_tmp2, Potential_tmp, thrust::plus<Real>());
thrust::transform(thrust::device, F_tmp, F_tmp+N-2, F_tmp2, F_tmp, thrust::plus<vec3>());

/*    
    vec3 FFF;
    Real PPP;
    SCFObject1.CalculateGravityAtPoint(vec3(100,0.001,0.001), &PPP, &FFF);
    
    std::cerr << PPP << std::endl;
    std::cerr << FFF.x << std::endl;
    std::cerr << FFF.y << std::endl;
    std::cerr << FFF.z << std::endl;*/

//     vec3 FFF;
//     Real PPP;
//     cudaMemcpy(&PPP, Potential_tmp+N-2, 1*sizeof(Real), cudaMemcpyDeviceToHost);
//     cudaMemcpy(&FFF, F_tmp+N-2, 1*sizeof(vec3), cudaMemcpyDeviceToHost);
//     std::cerr << "===" << std::endl;
//     std::cerr << PPP << std::endl;
//     std::cerr << FFF.x << std::endl;
//     std::cerr << FFF.y << std::endl;
//     std::cerr << FFF.z << std::endl;
//     cudaMemcpy(&PPP, Potential_tmp+N-1, 1*sizeof(Real), cudaMemcpyDeviceToHost);
//     cudaMemcpy(&FFF, F_tmp+N-1, 1*sizeof(vec3), cudaMemcpyDeviceToHost);
//     std::cerr << "+++" << std::endl;
//     std::cerr << PPP << std::endl;
//     std::cerr << FFF.x << std::endl;
//     std::cerr << FFF.y << std::endl;
//     std::cerr << FFF.z << std::endl;

    cudaMemcpy(Potential, Potential_tmp, N*sizeof(Real), cudaMemcpyDeviceToDevice);
    cudaMemcpy(F,         F_tmp,         N*sizeof(vec3), cudaMemcpyDeviceToDevice);
}