// Interface for LAPACK function

# ifndef LAPACK_INTER_H
# define LAPACK_INTER_H

#include <complex>
typedef int lapack_int;
typedef complex<float> lapack_complex_float;
typedef complex<double> lapack_complex_double;

#ifdef USE_LAPACK

  #define DGETRF dgetrf_
  #define DGETRS dgetrs_
  #define DGETRI dgetri_
  #define ZPPTRF zpptrf_
  #define ZPPTRI zpptri_
  #define ZHPEVD zhpevd_

  #ifdef __cplusplus
  extern "C" {
  #endif /* __cplusplus */
  void dgetrf_( const lapack_int* m, const lapack_int* n, double* a, const lapack_int* lda, 
                lapack_int* ipiv, lapack_int* info );
  void dgetrs_( const char* trans, const lapack_int* n, const lapack_int* nrhs, 
							  const double* a, const lapack_int* lda, const lapack_int* ipiv, 
							  double* b, const lapack_int* ldb, lapack_int* info );
  void dgetri_( const lapack_int* n, double* a, const lapack_int* lda, 
                const lapack_int* ipiv, double* work, const lapack_int* lwork, 
                lapack_int* info );
  void zpptrf_( const char* uplo, const lapack_int* n, lapack_complex_double* ap, 
                lapack_int* info );
  void zpptri_( const char* uplo, const lapack_int* n, lapack_complex_double* ap, 
                lapack_int* info );
  void zhpevd_( const char* jobz, const char* uplo, const lapack_int* n,
							  lapack_complex_double* ap, double* w, lapack_complex_double* z,
							  const lapack_int* ldz, lapack_complex_double* work, const lapack_int* lwork,
							  double* rwork, const lapack_int* lrwork, lapack_int* iwork,
							  const lapack_int* liwork, lapack_int* info );
  #ifdef __cplusplus
  }
  #endif /* __cplusplus */
  
#else

  //#define MKL_Complex8 lapack_complex_float
  //#define MKL_Complex16 lapack_complex_double
  #include "mkl.h"

  inline void ZPPTRF( char* uplo, const lapack_int* n, lapack_complex_double* ap, lapack_int* info ) {
    ZPPTRF(uplo, (int*)n, (MKL_Complex16*)ap, (int*)info);
  }
  inline void ZPPTRI( char* uplo, const lapack_int* n, lapack_complex_double* ap, lapack_int* info ){
     ZPPTRI(uplo, (int*)n, (MKL_Complex16*)ap, (int*)info);
  }

#endif

#endif /* lapack_intER_H */
