# ifndef __LINSYS_H
# define __LINSYS_H
typedef float *Pfloat;

# ifdef __cplusplus
extern "C" {
# endif

/* vectors are in columns of matr*/

typedef double lstype;
typedef lstype *lstypeP;      

int linsys(lstypeP *matr,lstype *vect,int num,lstype *x);

# ifdef __cplusplus
}
# endif

# endif

