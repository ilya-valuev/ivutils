// Copyright 2014, Nikita Kazeev, JIHT RAS
#include "../include/diff.h"

using std::fabs;
const double h = 1e-6;
const double diff_error_margin = 3e-2; // Used in complex diff existence verification


// TODO(kazeevn) switch to templates.
// I tried, but had problems with type resolving, as we need a separate function for
// complex diff

template<typename F>
inline double sdiff(const F& func, const double x) {
  return (-func(x+2*h)+8.*func(x+h)-8.*func(x-h)+func(x-2*h))/12./h;
}

inline double gsl_adaptor(const double x, void* func) {
  const auto F = *(static_cast<const std::function<double(double)>*>(func));
  return F(x);
}

double diff(std::function<double(double)> func, const double x, double* abserr) {
  double result;
  const gsl_function F = {&gsl_adaptor, static_cast<void*>(&func)};
  const int status = gsl_deriv_central(&F, x, h, &result, abserr);
  if (status != 0) {
    printf ("error: %s\n", gsl_strerror(status));
  }
  assert(status == 0);
  if (std::isnan(result))
    std::cerr << "\nx: " << x << "\nDiff: " << sdiff(func, x) << '\n';
  assert(!std::isnan(result));
  return result;
}

cdouble diff(const std::function<cdouble(double)>& func, const double x, double* abserr) {
  // There might be different optimal step sizes for real and diff
  // parts. Thus we calculate them separately.
  std::function<double(double)> func_real([&func] (const double xx) {return real(func(xx));});
  std::function<double(double)> func_imag([&func] (const double xx) {return imag(func(xx));});
  const gsl_function F_real = {&gsl_adaptor, static_cast<void*>(&func_real)};
  const gsl_function F_imag = {&gsl_adaptor, static_cast<void*>(&func_imag)};
  double result_real, abserr_real;
  const int status_real = gsl_deriv_central(&F_real, x, h, &result_real, &abserr_real);
  double result_imag, abserr_imag;
  const int status_imag = gsl_deriv_central(&F_imag, x, h, &result_imag, &abserr_imag);
  if (status_real != 0)
    printf ("error: %s\n", gsl_strerror(status_real));
  if (status_imag != 0)
    printf ("error: %s\n", gsl_strerror(status_imag));
  assert(status_real == 0);
  assert(status_imag == 0);
  assert(!std::isnan(result_real));
  assert(!std::isnan(result_imag));
  *abserr = abserr_real + abserr_imag;
  // std:cerr << "D :" << result_imag << " E: " << abserr_imag << " O: " <<
  //      (-func_imag(x+2*h)+8.*func_imag(x+h)-8.*func_imag(x-h)+func_imag(x-2*h))/12./h << '\n';
  return cdouble(result_real, result_imag);
}

cdouble diff(const std::function<cdouble(cdouble)>& func, const cdouble& x, double* abserr) {
  double abserr_real;
  const cdouble real_diff = diff(std::function<cdouble(double)>(
      [&func, &x] (const double xx) {
        return func(cdouble(xx, imag(x)));}), real(x), &abserr_real);
  double abserr_imag;
  const cdouble imag_diff = diff(std::function<cdouble(double)>(
      [&func, &x] (const double xx) {
        return func(cdouble(real(x), xx));}), imag(x), &abserr_imag);

  // Test, whether the derivative exists, using the Cauchy-Riemann equation
  // std::cerr << '\n' << real_diff << " pm " << abserr_real << '\n' <<
  //     imag_diff << " pm " << abserr_imag << '\n';
  const unsigned char TOLERANCE = 3;
  assert(std::fabs(real_diff.real() - imag_diff.imag()) < TOLERANCE * (
         real(abserr_real) + imag(abserr_imag) +
         diff_error_margin*(fabs(real_diff.real()) + fabs(imag_diff.imag()))));
  assert(std::fabs(real_diff.imag() + imag_diff.real()) < TOLERANCE * (
         imag(abserr_real) + real(abserr_imag) +
         diff_error_margin*(fabs(real_diff.imag()) + fabs(imag_diff.real()))));
  *abserr = abserr_real;
  return real_diff;
}
