/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author  : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project  : GridMD, ivutils
 *
 *
 *
 *****************************************************************************/
/*r @file vector_3.h @brief работа с трехмерными комплексными векторами
*/

# ifndef CVECTOR_3A_H
# define CVECTOR_3A_H

# include <complex>
# include <cmath>
# include "vector_3.h"
# include <iostream>

using namespace std;

typedef complex<vec_type>     cdouble;
typedef Vector_Nt<cdouble,3>  cVector_3;
//typedef Vector_Nt<cdouble, 0, ptr_stor_t<cdouble> > cVector_G;

//------------------------------------------------------
// Overloads for cdouble
//------------------------------------------------------

inline cdouble operator*(int a, const cdouble &b){
  return ((double)a)*b;
}
inline cdouble operator*(const cdouble &b,int a){
  return a*b;
}
inline cdouble operator/(const cdouble &b,int a){
  return (1./a)*b;
}

inline cdouble operator/(int a, const cdouble &b){
  return (a)*(1./b);
}

//------------------------------------------------------
// Overloads for cVector_3
//------------------------------------------------------
inline std::ostream& operator<<(::std::ostream& os, const cVector_3& vector) {
  return os << vector[0] << '\n' <<  vector[1] << '\n' << vector[2];
  }

inline cVector_3 operator*(const cdouble &a, const Vector_3 &v){
  return cVector_3(a*v[0], a*v[1], a*v[2]);  // a*cVector_3(v);
}

inline cVector_3 operator*(const Vector_3 &v, const cdouble &a){
  return a*v;
}

inline Vector_3 real(const cVector_3 &cv){
  return Vector_3(cv[0].real(), cv[1].real(), cv[2].real());
}

inline Vector_3 imag(const cVector_3 &cv){
  return Vector_3(cv[0].imag(), cv[1].imag(), cv[2].imag());
}

inline cVector_3 conj(const cVector_3 &cv){
  return cVector_3(conj(cv[0]), conj(cv[1]), conj(cv[2]));
}

inline cVector_3 rcell1(cVector_3& cv, Vector_3 &cell, int flags=0xffff) {
  return cVector_3( real(cv).rcell1(cell, flags) ) + cdouble(0,1)*imag(cv);
}

inline cVector_3 from_parts(const Vector_3 real_part, const Vector_3 imag_part) {
  cVector_3 res;
  for (unsigned char axis = 0; axis < 3; ++axis) {
    res[axis] = cdouble(real_part[axis], imag_part[axis]);
  }
  return res;
}

# endif // __CVECTOR_3A_H
