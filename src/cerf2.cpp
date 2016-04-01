/** @file cerf2.cpp 
    @brief Optimized complex error function. */
  

/*
$Source: /home/plasmacvs/source_tree/ivutils/src/cerf2.cpp,v $
$Revision: 1.1 $
$Author: morozov $
$Date: 2011/05/24 23:16:31 $
*/

//!!!DEBUG
//extern int cnt[];

#define XY_INTERPOLATION_TABLE

#include "cerf2.h"
namespace Octave {
#include "cerf_octave.h"
}

using namespace std;

const double CerfDataContainer::eps_min = 1e-13;
const double CerfDataContainer::eps_max = 1e-2;

CerfDataContainer cerf_data;  // global variable

int CerfDataContainer::init(double eps){
  int i;

  // Check eps range
  if(eps < eps_min || eps > eps_max) return 1;

  // ======= Tune cerf_continued_fraction ===========
  Octave::eps = pow(eps, 1.006880707) * 0.00011;
  Octave::eps2 = Octave::eps * Octave::eps;

  // ======= Localize the region where erf = 1 ===========
  // Calculate min/max rho for the table
  r0 = 1.02 * sqrt( -log(1.4*eps) ) - 0.285;
  if(r0 < 2.) r0 = 2.;  // compatibility with test abs(z) < 2. for series
  rr0 = r0*r0;
  double rr1 = 16.*rr0;   // r1 = 10*r0
  double rr_step = (rr1 - rr0) / (nrr - 2);
  rr_step_1 = 1. / rr_step;

  // Create the table phi(rho)
  xx_rr = new double[nrr];
  eps_pi = eps*eps*M_PI;
  for(i = 0; i < nrr; i++) {
    double rr = rr0 + i*rr_step;  // maximal rr for the interval [i;i+1]
    xx_rr[i] = (rr - log( sqrt(rr*eps_pi) )) / 2.;
  }
  cos2_phi1 = xx_rr[nrr1] / rr1;

#ifdef XY_INTERPOLATION_TABLE
  // ======= XY interpolation table for rho < rho0 ===========
  double delta = pow(eps, 0.465) * 6.4;
  if(delta > 0.02) delta = 0.02;
  
  x0 = 2.;
  double x1 = r0, y1 = 2.;
  double dx = delta, dy = delta;
  nx = (int)((x1 - x0) / dx) + 2;
  ny = (int)(y1 / dy) + 2;
  if((double)nx*ny > 8192*4096) return 2;  // Table is too large!

  if(nx > 1 && ny > 1) {
    tab = new cdouble[nx*ny]; 
    nx1 = nx - 1; ny1 = ny - 1;
    xstep_1 = 1. / dx, ystep_1 = 1. / dy;

    for(i = 0; i < nx; i++) {
      double x = x0 + i*dx;
      cdouble *tabx = tab + i*ny;
      for(int j = 0; j < ny; j++) {
        cdouble z(x, j*dy);
        tabx[j] = (1. - cerf_Octave_exact(z)) / (exp(-z*z)/z);
      }
    }
  }
  else nx = nx1 = ny = ny1 = 0;
#endif

  // ======= rho-phi interpolation table for rho < rho0 ===========
  /*nr = 2048; np = 2048;
  nr1 = nr - 1; np1 = np - 1;
  tab = new cdouble[nr*np]; 

  double dr = r0 / (nr - 2);  rstep_1 = 1. / dr;
  double dp = M_PI / 2. / (np - 2);  pstep_1 = 1. / dp;
  for(i = 0; i < nr; i++) {
    double rho = i*dr;
    cdouble *tabr = tab + i*np;
    for(int j = 0; j < np; j++) {
      double phi = j*dp;
      tabr[j] = cerf_Octave_exact( cdouble(rho*cos(phi), rho*sin(phi)) );
    }
  }*/

  return 0;
}

int cerf_data_init(double eps) { 
  return cerf_data.init(eps);
}


inline cdouble unity() {  //!!!DEBUG
  //cnt[1]++;
  return 1.;
}


cdouble cerf(const cdouble z) {
  cdouble res;

#ifdef XY_INTERPOLATION_TABLE
  // Bilinear interpolation using XY table
  cdouble za( fabs(z.real()), fabs(z.imag()) );
  double xt = (za.real() - cerf_data.x0) * cerf_data.xstep_1;
  double yt = za.imag() * cerf_data.ystep_1;
  int ix = (int)xt,  iy = (int)yt;
  if(ix >= 0 && ix < cerf_data.nx1 && iy < cerf_data.ny1) {
    xt -= ix;  yt -= iy;
    int ny = cerf_data.ny;
    int idx = ix*ny + iy;
    res = cerf_data.tab[idx]*((1-xt)*(1-yt)) + cerf_data.tab[idx+ny]*(xt*(1-yt)) +
          cerf_data.tab[idx+1]*((1-xt)*yt)   + cerf_data.tab[idx+ny+1]*(xt*yt);
    res = 1. - res * exp(-za*za) / za;

    if(z.real() < 0) res = - conj(res);
    if(z.imag() < 0) res = conj(res);
    return res;
  }
#endif

  // Bilinear interpolation using rho-phi table
  /*cdouble za( fabs(z.real()), fabs(z.imag()) );
  double r = abs(za) * cerf_data.rstep_1, p;
  int ir = (int)r, ip;
  if(ir < cerf_data.nr1) {
    r -= ir;
    p = arg(za) * cerf_data.pstep_1;
    ip = (int)p;
    p -= ip;
    int np = cerf_data.np;
    int idx = ir*np + ip;
    res = cerf_data.tab[idx]*((1-r)*(1-p)) + cerf_data.tab[idx+np]*(r*(1-p)) +
          cerf_data.tab[idx+1]*((1-r)*p)   + cerf_data.tab[idx+np+1]*(r*p);

    if(z.real() < 0) res = - conj(res);
    if(z.imag() < 0) res = conj(res);
    return res;
  }*/

  double x = z.real(), xx = x*x;
  double rr = xx + z.imag()*z.imag();

  if (rr < 4.)
    res = Octave::cerf_series( z );
  else if (fabs(x) < 0.5)
    res = Octave::cerf_rybicki( z );
  else {
    //cnt[0]++;
    // Process the case of erf ~ 1
    //if( xx > .5*( rr - log( sqrt(rr*cerf_data.eps_pi) ) ) )  res = unity();
    if(rr > cerf_data.rr0) {
      double drr = (rr - cerf_data.rr0) * cerf_data.rr_step_1;
      int irr = (int)drr;
      if(irr < cerf_data.nrr1) {
        drr -= irr;
        res = ( xx > (cerf_data.xx_rr[irr]*(1. - drr) + cerf_data.xx_rr[irr+1]*drr) ) ?
              unity() : Octave::cerf_continued_fraction( z );
      }
      else {
        //cnt[2]++;
        res = (xx > rr*cerf_data.cos2_phi1) ? unity() : Octave::cerf_continued_fraction( z );
      }
    }
    else res = Octave::cerf_continued_fraction( z );
  }
  
  return res;
}


cdouble cerf_Octave(const cdouble z){
  return Octave::cerf(z);
}


cdouble cerf_Octave_exact(const cdouble z){
  double eps = Octave::eps;
  Octave::eps = 1e-15, Octave::eps2 = Octave::eps * Octave::eps;
  cdouble res = Octave::cerf(z);
  Octave::eps = eps, Octave::eps2 = Octave::eps * Octave::eps;
  return res;
}
