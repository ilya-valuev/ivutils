/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.6 $
 *   $Date: 2014/04/17 13:51:00 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/vector_n.h,v 1.6 2014/04/17 13:51:00 kazeev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/vector_n.h,v $
$Revision: 1.6 $
$Author: kazeev $
$Date: 2014/04/17 13:51:00 $
*/
/*e****************************************************************************
 * $Log: vector_n.h,v $
 * Revision 1.6  2014/04/17 13:51:00  kazeev
 * Fixed gcc compilation. Added BoxHamiltonian.
 *
 * Revision 1.5  2009/03/04 09:52:45  valuev
 * corrected after sync with FDTD project
 *
 * Revision 1.4  2009/02/10 14:20:45  valuev
 * sync with FDTD project
 *
 * Revision 1.3  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.3  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.2  2006/08/08 13:02:04  valuev
 * Added geometry
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:11  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/
# if 0  //      ndef _Windows
# include <graphics.h>
# endif

# ifndef __VECTOR_N_H

# define __VECTOR_N_H


# include "common.h"

#define DOUBLEUSED_N

# ifdef DOUBLEUSED_N
 typedef double vec_typen;
# define PRTCHARN "d"
# else
 typedef float vec_typen;
# define PRTCHARN "f"
# endif

class Hermitian;

class Vector_Np {
 public:
  int size;
 private:
  void chksz(int sz){
    if(size!=sz)serror("Vector_Np: sizes do not match.\n");
  }
 public:
  int alloc;
  vec_typen *v;


  //e wrapper, copying the POINTER ONLY
  Vector_Np(vec_typen *vect, int sz){
    size=sz;
    v=vect;
    alloc=0;
  }

  //e copying the array if not NULL
  Vector_Np(int sz,vec_typen *vect=NULL){
    v=new vec_typen[sz];
    //v=(vec_typen *)malloc(sz*sizeof(vec_typen));
    if(!v)serror("Vector_Np: memory allocation error.\n");
    size=sz;
    if(vect){
      int i;
      for(i=0;i<size;i++)v[i]=vect[i];
    }

    alloc=1;
  }

  ~Vector_Np(){
    if(alloc)delete [] v;
  }

  inline vec_typen& operator[](int i){ return v[i]; };

  Vector_Np& operator=(float *vect);
  Vector_Np& operator=(double *vect);

  Vector_Np& operator=(Vector_Np& vect);

  Vector_Np operator+(Vector_Np& vect);
  Vector_Np operator-(Vector_Np& vect);


  vec_typen operator*(Vector_Np& vect);

  Vector_Np operator*(vec_typen coeff);
  Vector_Np operator/(vec_typen coeff);

  friend Vector_Np operator*(vec_typen coeff,Vector_Np& vec);

  Vector_Np operator-();

  Vector_Np& operator*=(vec_typen coeff);
  Vector_Np& operator/=(vec_typen coeff);

  Vector_Np& operator+=(Vector_Np& vect);
  Vector_Np& operator-=(Vector_Np& vect);


  vec_typen norm();
  vec_typen normalize();

  void print(char *form=NULL);
  operator vec_typen*(){
    return v ;
  }

  friend Vector_Np operator*(Hermitian &A, Vector_Np &vec);
  friend int BandMul(Vector_Np &resvec,int ns,vec_typen *band_mat,Vector_Np &vec);

};


// multiplies the band matrix of with n1-1 subdiagonals by vec
int BandMul(Vector_Np &resvec,int ns,vec_typen *band_mat,Vector_Np &vec);



// neue Funktionen



# endif // __VECTOR_N_H
