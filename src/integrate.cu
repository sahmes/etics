// this was a dumb idea, move everything back to main, and include main.hpp or something in the amuse interface

// only global variables, functions and structures related to the leapfrog integration; the force calculation routines are in mex.cu and scf.cu
// this is basically auxiliary to main.cu
#include <thrust/device_vector.h>
#include "common.hpp"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#warning Those variables should actually be in main.cu, not here
Particle *P_h;
thrust::device_vector<Particle> P;

thrust::device_vector<Real> Potential;
thrust::device_vector<vec3> F0;
thrust::device_vector<vec3> F1;

Real ConstantStep = 0.001953125;
Real T, Step, eta, dT1, dT2, Tcrit;
int N;


int CommitParticles() {
    P = thrust::device_vector<Particle>(P_h, P_h+N);
    return 0;
}

int InitilizeIntegratorMemory() {
    Potential = thrust::device_vector<Real>(N);
    F0 = thrust::device_vector<vec3>(N);
    F1 = thrust::device_vector<vec3>(N);
    return 0;
}

struct CommitForcesFunctor {
    Real Step;
    __host__ __device__ CommitForcesFunctor(Real _Step) : Step(_Step) {}
    __host__ __device__ Particle operator() (Particle& p, const vec3& F) const {
        p.acc = F;
        return p;
    }
};

void CommitForces() {
    thrust::transform(P.begin(), P.end(), F0.begin(), P.begin(), CommitForcesFunctor(Step));
}

// The 'drift' step is performed using the 'acc' member.
struct DriftFunctor {
    Real Step, Step2;
    __host__ __device__ DriftFunctor(Real _Step) : Step(_Step), Step2(_Step*_Step) {}
    __host__ __device__ Particle operator() (Particle &p) const {
        p.pos += p.vel*Step + p.acc*0.5*Step2;
        p.CalculateR2(); // needed for both MEX and SCF, but not generally needed in leapfrog
        return p;
    }
};

void DriftStep() {
    thrust::transform(P.begin(), P.end(), P.begin(), DriftFunctor(Step));
}

// The 'kick' step is performed using the 'acc' member and also the force F,
// calculated at the new (predicted) position.
struct KickFunctor {
    Real Step;
    __host__ __device__ KickFunctor(Real _Step) : Step(_Step) {}
    __host__ __device__ Particle operator() (Particle& p, const vec3& F) const {
        p.vel += (p.acc + F)*0.5*Step;
        p.acc = F;
        return p;
    }
};

void KickStep() {
    thrust::transform(P.begin(), P.end(), F1.begin(), P.begin(), KickFunctor(Step));
}
