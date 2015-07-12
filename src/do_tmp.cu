#include <iostream>
#define Real double
__host__ __device__ Real Pl(int l, Real x);

template <int l>
__host__ __device__ Real Pl_new(Real x);






template <int N> 
struct Factorial {
  static const int value = N * Factorial<N - 1>::value;
};
 
// Base case via template specialization:
 
template <>
struct Factorial<0> {
  static const int value = 1;
};






template <int N> 
struct Yohai {
  static const int value = 88;
};
 
// Base case via template specialization:


template <>
struct Yohai<0> {
  static const int value = 99;
};

template <>
struct Yohai<1> {
  static const int value = 99;
};

template <>
struct Yohai<2> {
  static const int value = 99;
};


template <>
struct Yohai<3> {
  static const int value = 99;
};


int main() {
//     double z = 0;
//     z += Factorial<0>::value;
//     z += Factorial<1>::value;
//     z += Factorial<2>::value;
//     z += Factorial<3>::value;
//     std::cout << z << std::endl;
    
    for (int l=0; l < 3; l++) {
        std::cout << Pl(l, 0.5) << std::endl;
        std::cout << Yohai<l>::value << std::endl;
//         std::cout << Factorial<6>::value << std::endl;
//         std::cout << Pl_new<l>(0.5) << std::endl;
    }
    
    
}