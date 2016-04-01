# ifndef __CVECTOR_3_H
# define __CVECTOR_3_H

# include <complex>

# include "vector_3.h"

using namespace std;

typedef complex<vec_type> cvec_type;

class cVector_3 {
  Vector_3 V[2];

public:
  cVector_3(const Vector_3 &re=Vector_3(), const Vector_3 &im=Vector_3()) {
    V[0]=re;
    V[1]=im;
  }

  inline Vector_3& re() { 
    return this->V[0];
  };

  inline Vector_3& im() { 
    return this->V[1];
  };

  inline Vector_3& real() { 
    return this->V[0];
  };

  inline Vector_3& imag() { 
    return this->V[1];
  };

  inline int operator==(const cVector_3 &cvect) const{
    return (V[0]==cvect.V[0] && V[1]==cvect.V[1]);
  };

  inline int operator!=(const cVector_3 &cvect) const{
    return (!(*this==cvect));
  };

  inline cVector_3 operator+(const cVector_3& cvect) const{
    return cVector_3(V[0]+cvect.V[0], V[1]+cvect.V[1]);
  }

  inline cVector_3 operator-(const cVector_3 &cvect) const {
    return cVector_3(V[0]-cvect.V[0], V[1]-cvect.V[1]);
  }

  inline cVector_3& operator+=(cVector_3 cvect){
    V[0]+=cvect.V[0];
    V[1]+=cvect.V[1];
    return *this;
  }

  inline cvec_type operator*(const cVector_3& cvect) const {
    return cvec_type(V[0]*cvect.V[0]-V[1]*cvect.V[1], V[0]*cvect.V[1]+V[1]*cvect.V[0]);
  }

  inline cVector_3 operator %(const cVector_3 &cvect) const{
    return cVector_3(V[0]%cvect.V[0]-V[1]%cvect.V[1], V[0]%cvect.V[1]+V[1]%cvect.V[0]);
  }

  inline cVector_3 operator %(const Vector_3 &vect) const{
    return cVector_3(V[0]%vect, V[1]%vect);
  }

  inline cVector_3 operator*(vec_type coeff) const {
    return cVector_3(V[0]*coeff, V[1]*coeff);
  }

  inline cVector_3 operator*(cvec_type coeff) const {
    return cVector_3(V[0]*coeff.real() - V[1]*coeff.imag(), V[1]*coeff.real() + V[0]*coeff.imag());
  }

  friend cVector_3 operator*(vec_type coeff,const cVector_3& cvec);
  friend cVector_3 operator*(cvec_type coeff,const cVector_3& cvec);

  inline cVector_3 operator/(vec_type coeff){
    return cVector_3(V[0]/coeff, V[1]/coeff);
  }

  inline cVector_3 operator-(){
    return cVector_3(-V[0], -V[1]);
  }

  inline cVector_3& operator*=(vec_type coeff){
    V[0]*=coeff;
    V[1]*=coeff;
    return *this;
  }
};

# endif
