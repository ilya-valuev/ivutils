// Copyright 2014, Nikita Kazeev, JIHT RAS
#ifndef DIFF_H
#define DIFF_H

//#include <type_traits>
#include <functional>
#include <assert.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include "./cvector_3.h"

double diff(std::function<double(double)> func, const double x, double* abserr);
cdouble diff(const std::function<cdouble(double)>& func, const double x, double* abserr);
cdouble diff(const std::function<cdouble(cdouble)>& func, const cdouble& x, double* abserr);

template<typename T_arg, typename T_val, unsigned int N>
Vector_Nt<T_val, N> diff_by_vector(
    std::function<T_val(Vector_Nt<T_arg, N>)> func,
    const Vector_Nt<T_arg,3>& x, double* abserr) {
  Vector_Nt<T_val, N> res;
  *abserr = 0;
  for (unsigned char axis = 0; axis < N; ++axis) {
    double axis_err;
    res[axis] = diff(std::function<T_val(T_arg)>(
        [&func, &x, axis] (const T_arg xx) -> T_val {
          Vector_Nt<T_arg, N> y(x);
          y[axis] = xx;
          return func(y);
        }), x[axis], &axis_err);
    *abserr += std::fabs(axis_err);
  }
  return res;
}


#endif // DIFF_H
