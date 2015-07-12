#pragma once
#include <string>

// A_ON_SHARED_MEMORY moves the A structure (after it has been calculated) from constant to shared memory. This might be faster, or not.
// #define A_ON_SHARED_MEMORY


#if defined SINGLE_PRECISION && defined DOUBLE_PRECISION
    #error Contradictory precision flags!
#endif

#ifndef SINGLE_PRECISION
    #define DOUBLE_PRECISION
    #define Real double
#else
    #define Real float
#endif

#if defined MEX && defined SCF
    #error Contradictory method flags!
#endif

#if defined MEX && !defined LMAX
    #error LMAX not defined.
#endif

#if defined SCF && ((!defined NMAX) || (!defined LMAX))
    #error Both NMAX and LMAX should be defined!
#endif

#if defined MEX && defined NMAX
    #warning NMAX is defined, but will be ignored since the method is MEX!
#endif

struct vec3 {
    Real x, y, z;

    __host__ __device__ vec3() : x(0), y(0), z(0) {}
    __host__ __device__ vec3(Real _x, Real _y, Real _z) : x(_x), y(_y), z(_z) {}
    __host__ __device__ Real abs2() const {return x*x + y*y + z*z;}
    __host__ __device__ vec3 operator+ (const vec3& V)   const {return vec3(this->x + V.x, this->y + V.y, this->z + V.z);}
    __host__ __device__ vec3 operator- (const vec3& V)   const {return vec3(this->x - V.x, this->y - V.y, this->z - V.z);}
    __host__ __device__ vec3 operator* (const Real& C) const {return vec3(this->x*C, this->y*C, this->z*C);}
    __host__ __device__ vec3& operator+= (const vec3& V) {this->x += V.x; this->y += V.y; this->z += V.z; return *this;}
    __host__ __device__ vec3& operator-= (const vec3& V) {this->x -= V.x; this->y -= V.y; this->z -= V.z; return *this;}
    __host__ __device__ vec3 operator- () const {return vec3(-this->x, -this->y, -this->z);}
    __host__ __device__ vec3 operator+ () const {return vec3(this->x, this->y, this->z);}
};


class Particle {
    public:
        int ID;
        unsigned char Status;
        Real m;
        vec3 pos, vel, acc;
        Real R2;

        __host__ __device__ Particle() {}
        __host__ __device__ void CalculateR2() {R2 = pos.abs2();}
        __host__ __device__ bool operator< (const Particle& p) const {return (this->R2 < p.R2);}
};

struct ParametersStruct {
    int N;
    std::string FileName;
    int Skip;
    Real dT1, dT2, Tcrit, ConstantStep;
    std::string Prefix;
    int DeviceID;
};
