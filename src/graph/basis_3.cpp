#include "basis_3.h"
#include "../src/linsolv.hpp"

template<int N>
typename Basis_N<N>::vector_t Basis_N<N>::operator()(const vector_t &x) const{
  vector_t y;
  for(int i=0;i<N;i++){
    y[i]=0;
    for(int j=0;j<N;j++){
      y[i]+=vec[j][i]*x[j];
    }
  }
  return y;
}

template<int N>
typename Basis_N<N>::vector_t Basis_N<N>::cov(const vector_t &x) const{
  vector_t y;
  for(int i=0;i<N;i++){
    y[i]=0;
    for(int j=0;j<N;j++){
      y[i]+=vec[i][j]*x[j];
    }
  }
  return y;
}

template<int N>
typename Basis_N<N>::vector_t Basis_N<N>::inv(const vector_t &x) const{
  vec_type matr[N*N];
//  vec_type *pm[N];
  vector_t y, xx;
  for(int i=0;i<N;i++){
//    pm[i]=matr+i*N;
    xx[i]=x[i];
    for(int j=0;j<N;j++)
      matr[i*N+j]=vec[i][j];
//      pm[i][j]=vec[i][j];//(*this)(j,i);
  }
//  linsys_cpp(pm,xx.v,N,y.v);
//  linsys_cpp(matr,xx.get_ptr(),N,y.get_ptr());
    linsys(matr,xx.get_ptr(),N,y.get_ptr());
  //Vector_3 chk=(*this)(y);
  return y;
}

template<int N>
Basis_N<N> Basis_N<N>::inv() const{
  Basis_N b;
  for(int i=0;i<N;i++){
    Vector_3 x;
    x[i]=1;
    b.vec[i]=inv(x);
  }
  return b;
}

template class Basis_N<3>;

//template int linsys_cpp<double>(double **,double*, int, double *, int);