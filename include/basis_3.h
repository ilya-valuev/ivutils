/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 1995-2009        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.5 $
 *   $Date: 2013/08/16 11:32:45 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/basis_3.h,v 1.5 2013/08/16 11:32:45 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/basis_3.h,v $
$Revision: 1.5 $
$Author: valuev $
$Date: 2013/08/16 11:32:45 $
*/
/*s****************************************************************************
 * $Log: basis_3.h,v $
 * Revision 1.5  2013/08/16 11:32:45  valuev
 * compiled tcpengine
 *
 * Revision 1.4  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.10  2012/03/29 11:35:09  valuev
 * universal scatter_data() for emSourceWave
 *
 * Revision 1.9  2012/03/21 17:09:37  lesha
 * documentation
 *
 * Revision 1.8  2011/03/04 16:53:04  valuev
 * merged with gpu branch
 *
 * Revision 1.7  2011/01/26 00:59:12  lesha
 * TestLine is added to Inverse and StretchedRegion
 *
 * Revision 1.6  2010/10/21 16:09:10  biaks
 * add brief
 *
 * Revision 1.5  2010/01/02 22:22:37  lesha
 * StretchedRegion is added, Contour::GetCenter if fixed, 2DDumping is added etc.
 *
 * Revision 1.4  2010/01/01 20:40:07  lesha
 * Region, Box_N, Sphere_N, Basis_N, MonteCarloTestContour are added
 *
 * Revision 1.3  2009/09/30 14:44:57  valuev
 * Added fill_lattice function
 *
 * Revision 1.2  2009/06/04 14:33:42  lesha
 * *** empty log message ***
 *
 * Revision 1.1  2009/06/01 13:03:53  valuev
 * separated Basis_3 from Vector_3
 *
*******************************************************************************/
/// \en @file vector_3.h \brief 3D Basis
/// \ru @file vector_3.h \brief 3D базис

#ifndef BASIS_3_H

#define BASIS_3_H

#include "vector_3.h"

template<int N>
class Basis_N{
public:
  typedef Vector_Nt<vec_type,N> vector_t;
  vector_t vec[N];
  /* basis vectors are in columns, as in linear algebra */
/* Basis_3(vec_type b00,vec_type b01=0,vec_type b02=0,
         vec_type b10=0,vec_type b11=0,vec_type b12=0,
         vec_type b20=0,vec_type b21=0,vec_type b22=0){
  vec[0][0]=b00;vec[1][0]=b01;vec[2][0]=b02;
  vec[0][1]=b10;vec[1][1]=b11;vec[2][1]=b12;
  vec[0][2]=b20;vec[1][2]=b21;vec[2][2]=b22;
 }*/
  Basis_N(const vector_t &x1, const vector_t &x2, const vector_t &x3){
    if(N>0)vec[0]=x1;
    if(N>1)vec[1]=x2;
    if(N>2)vec[2]=x3;
  }
  Basis_N(vec_type stretch=1.){
    for(int i=0;i<N;i++)
      vec[i][i]=stretch;
  }
  Basis_N(const Basis_N &other){
    for(int i=0;i<N;i++)
      vec[i]=other.vec[i];
  }
  
  vec_type &operator()(int i, int j){
    return vec[j][i]; 
  }

  vec_type operator()(int i, int j) const{
    return vec[j][i]; 
  }

 /// transforms a vector from the orthogonal basis to this basis
 /// (contravariant form y=B*x)
  vector_t operator()(const vector_t &x) const;

 /// transforms a vector from the orthogonal basis to this basis
 /// (covariant form y=transpose(B)*x )
  vector_t cov(const vector_t &x) const;
  
 /// transforms a vector from this basis into the orthogonal basis
 /// (contravariant, x=B^(-1)*y )
  vector_t inv(const vector_t &x) const;

  /// returns inverse basis
  Basis_N inv() const;

  inline vector_t& operator[](int i) {
    return vec[i];
  }

  inline vector_t& operator[](size_t i) {
    return vec[i];
  }

  inline vector_t operator[](int i) const {
    return vec[i];
  }

  inline vector_t operator[](size_t i) const {
    return vec[i];
  }

  inline vec_type volume() const {
    return (vec[0]%vec[1]%vec[2]).norm();
  }
};

typedef Basis_N<3> Basis_3;

template<int N>
Basis_N<N> get_stretch_basis(vec_type stretch){
  return Basis_N<N>(stretch);	
  /*Basis_N<N> b;
  for(int i=0;i<N;i++)
    b.vec[i][i]=stretch;
  return b;*/
}


/// basic class for vector transformation (default is no change)
struct VecTransform{
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return arg;
  }
  virtual VecTransform *copy() const {
    return new VecTransform();
  }
  virtual ~VecTransform(){}
};

/// vector constant shift
struct VecShift: public VecTransform{
  Vector_3 shift;
  VecShift(const Vector_3 &shift_=Vector_3()):shift(shift_){}
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return arg+shift;
  }
  virtual VecTransform *copy() const {
    return new VecShift(shift);
  }
};

/// vector homothety
struct VecHomothety: public VecTransform{
  vec_type a;
  VecHomothety(const vec_type a_):a(a_){}
  virtual Vector_3 operator()(const Vector_3& arg) const{
    return a*arg;
  }
  virtual VecTransform *copy() const {
    return new VecHomothety(a);
  }
};

#endif
