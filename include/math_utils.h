#ifndef _MATH_UTILS_H
#define _MATH_UTILS_H

/// \en @file math_utils.h \brief Some useful math functions

#include <functional>
#include <complex>
#include "refobj.h"
#include "vector_3.h"

/// tests if x is NaN
template<class T>
int test_NaN(T x){
  return (!(x>0 || x<0 || x==0));
}

///\en round given value up to the near number of the type 2^N
///\ru округляет число вверх до ближайшего числа вида 2^N
template<class T>
T ceil2(const T val){
  int iexp;
  frexp(val,&iexp);
  T res=ldexp(.5, iexp+1);
  return res;
}

///\en compare x and y with some accuracy which depends on parameter val.
///\ru сравнивает x и y с точностью до некоторого малого числа (зависящего от параметра val)
///    по умолчанию y=0
template<class T>
bool acccomp(const T x, const T y, const T val=T(1), int mult_epsilon=1024) {
  T mult=ceil2(val);
  T err=mult*T(mult_epsilon)*numeric_limits<T>::epsilon();
  return fabs(x-y)<=err;
}

template<class T>
bool accless(const T x, const T y, const T val=T(1)) {
  if(x<=y || acccomp(x,y,val))return true;
  return false;
}

///\en calculates ratio between x and y.
/// if this ratio is close to some integer number with some accuracy,
/// then returns this integer number.
///\ru делит x на y. Если их частное меньше некоторого малого числа (зависящего от величины частного), то оно округляется.
template<class T>
T accdiv(T x, T y, int mult_epsilon=1024) {
  T fr=x/y;
  T ifr=floor(fr+T(.5));
  int mult = fabs(ifr)<=T(mult_epsilon) ? mult_epsilon : int(ceil2(fabs(ifr)));
  T err=mult*numeric_limits<T>::epsilon();
  if (fabs(fr-ifr)<=err)
    return ifr;
  else
    return fr;
}

///\en Logarithmically reduces the argument:
/// log_reduce(x) = log(x), x>e
/// log_reduce(x) -> -log(-x), x<e
/// linear between [-e,e], first derivative continuous
template <class T>
T log_reduce(const T& x){
  return x<-M_E ? -log(-x) : x<M_E ? x/M_E : log(x) ;
  //return x<0 ? -x*x : x*x ;
  //return x<0 ? -sqrt(-x) : x==0 ? 0 : sqrt(x) ;
  //return x;
}

template<class T>
T real_value(const T &a){
  return a;
}

template<class T>
T real_value(const complex<T> &a){
  return a.real();
}

/// returns -1 if t1<t2; 1 if t2>t1; 0 if t1==t2
/// can be used in implementation of operators < or > for some classes
template<class T>
int step_fun(const T &t1, const T &t2){
  if(t1<t2)return -1;
  else if(t2<t1)return 1;
  return 0;
}

template<class T>
int step_fun(const complex<T> &a1, const complex<T> &a2){
  if(a1.real()<a2.real())return -1;
  else if(a1.real()>a2.real())return 1;
  else if(a1.imag()<a2.imag())return -1;
  else if(a1.imag()>a2.imag())return 1;
  return 0;
}

template<class value_t, int N>
int step_fun(const Vector_Nt<value_t, N> &t1, const Vector_Nt<value_t, N> &t2){
  for(int i=0;i<N;i++){
    int res=step_fun(t1[i],t2[i]);
    if(res)return res;
  }
  return 0;
}

/*
template<template<class T> class containter_r>
int step_fun(const containter_r<class T> &t1, const containter_r<class T> &t2){
  if(t1.size()<t2.size())return -1;
  if(t1.size()>t2.size())return 1;

  for(typename containter_r::const_iterator it1=t1.begin(),it2=t2.begin(),e=t1.end();it1!=e;++it1,++it2){
    int res=step_fun(*it1,*it2);
    if(res)return res;
  }
  return 0;
}
*/
/*
template<class T>
int step_fun(const vector<T> &t1, const vector<T> &t2){
  if(t1.size()<t2.size())return -1;
  if(t1.size()>t2.size())return 1;

  for(typename vector<T>::const_iterator it1=t1.begin(),it2=t2.begin(),e=t1.end();it1!=e;++it1,++it2){
    int res=step_fun(*it1,*it2);
    if(res)return res;
  }
  return 0;
}
*/

template<class arg_t, class result_t>
class virt_unary_function: public unary_function<arg_t,result_t>{
public:
  virtual result_t operator()(arg_t x){return 0;}
};

template<class arg_t, class result_t>
class delta_function: public virt_unary_function<arg_t,result_t>{
  arg_t x0;
public:
  delta_function(arg_t x0_):x0(x0_){}
  virtual result_t operator()(arg_t x){
    return acccomp(x,x0) ? 1 : 0;
  }
};

/// used to construct max_function and min_function
template<class arg_t, class result_t, class comp_t=greater<result_t> >
class comp_function: public virt_unary_function<arg_t,result_t>{
protected:
  typedef virt_unary_function<arg_t,result_t> fun_t;
  refvector<fun_t> funs;
public:
  void add_function(fun_t *fun){
    funs.push_back(fun);
  }
  virtual result_t operator()(arg_t x){
    int sign = comp_t()(1,-1) ? 1 : -1;
    result_t res=-sign*numeric_limits<result_t>::max();
    for(size_t i=0;i<funs.size();i++){
      result_t y=funs[i]->operator()(x);
      comp_t()(y,res);
      if(sign>0 ? y>res : y<res)
        res=y;
    }
    return res;
  }
};

/// maximal value of given functions
template<class arg_t, class result_t>
class max_function: public comp_function<arg_t,result_t,greater<result_t> >{};

/// minimal value of given functions
template<class arg_t, class result_t>
class min_function: public comp_function<arg_t,result_t,less<result_t> >{};


///\en Pointer storage class for a vector followinfg \ref storage_prototype.
template <class T>
class ptr_stor_t{
  ///\en The data starts at a location in memory with specified address.
  shptr<T, delete_arr<T> > data;
  size_t dim_;
public:
  typedef T value_type;

  ///\en Default constructor.
  ptr_stor_t(): dim_(0){}

  ///\en Copy constructor.
  ptr_stor_t(const ptr_stor_t & other): data(other.data), dim_(other.dim_){}

  ///\en Binds the address in memory, if \a managed is 1, ownership over the pointer is transferred to the object
  ///    (the pointer will be deleted automatically).
  ptr_stor_t(size_t dimension, T *ptr=NULL, int managed=0): dim_(dimension){
    if(dim_){
      if(ptr)
        data.reset(ptr, managed);
      else
        data.reset(new T[dim_],1);
    }
  }

  ///\en Binds to existing std::vector
  ptr_stor_t(const std::vector<T> &vec): dim_(vec.size()){
    if(dim_)
      data.reset((T *)&(vec[0]),0);
  }

  ///\en Resize creates a new pointer, discarding the old one regardless of its size.
  void resize(size_t new_dim){
    data.reset(new T[new_dim],1);
    dim_=new_dim;
  }

  void resize(size_t new_dim, const T& value){
    resize(new_dim);
    for(size_t i=0;i<new_dim;i++)
      data[i]=value;
  }

  inline size_t dim() const { return dim_; }

  inline T& operator[](int i) const {return data.ptr()[i];}

  inline T& operator[](size_t i) const {return data.ptr()[i];}
};

template <class T>
ptr_stor_t<T> make_ptr_stor(const std::vector<T> &vec){
  return ptr_stor_t<T>(vec);
}

template <class T>
ptr_stor_t<T> make_ptr_stor(size_t dimension, T *ptr=NULL, int managed=0){
  return ptr_stor_t<T>(dimension, NULL, 0);
}

///\en General variable dimension vector around a pointer, N is redundant
typedef Vector_Nt<vec_type, 0, ptr_stor_t<vec_type> > Vector_G;


#endif
