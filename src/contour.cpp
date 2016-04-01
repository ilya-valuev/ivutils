/// \file \brief Non-template function definitions for  contour.h

#include "contour.h"
#include "pencil.h"

int set_perpendiculars(const Vector_3 &z, Vector_3 &x, Vector_3 &y){
  int i=0;
  for(;i<3;i++){
    if(z[i])break;
  }
  if(i==3)return -1;
  int j=(i+1)%3;
  x[i]=z[j];
  x[j]=-z[i];
  x[(j+1)%3]=0;
  x.normalize();
  y=z%x;
  return 0;
}

pencil<Vector_3> _PStmpv;
pencil<int> _PStmps;
