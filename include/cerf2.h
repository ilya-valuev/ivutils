# ifndef CERF2_H
# define CERF2_H

/** @file cerf2.h 
    @brief Header for the optimized complex error function. */
  

/*
$Source: /home/plasmacvs/source_tree/ivutils/include/cerf2.h,v $
$Revision: 1.1 $
$Author: morozov $
$Date: 2011/05/24 23:16:31 $
*/

#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>

using namespace std;

//#include "cvector_3.h"
typedef complex<double> cdouble;

int cerf_data_init(double eps);
cdouble cerf(const cdouble z);
cdouble cerf_Octave(const cdouble z);
cdouble cerf_Octave_exact(const cdouble z);


class CerfDataContainer {
public:
  static const double eps_min, eps_max;

  double *xx_rr;
  static const int nrr = 512, nrr1 = nrr - 1;
  double r0, rr0, rr_step_1, cos2_phi1, eps_pi;

  cdouble *tab;
  int nx, nx1, ny, ny1;
  double x0, xstep_1, ystep_1;
  /*int nr, nr1, np, np1;
  double rstep_1, pstep_1;*/

  CerfDataContainer() : xx_rr(NULL), tab(NULL) {}
  ~CerfDataContainer() {
    if(xx_rr) delete[] xx_rr;
    if(tab) delete[] tab;
  }
  int init(double eps);
};

extern CerfDataContainer cerf_data;  // global variable

# endif
