#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include "region_2.hpp"

Polygon_2 *GetRegularPolygon(const Vector_2 &center, vec_type d, int n){
  VecContour<2> *cnt = new VecContour<2>();
  vec_type dfi=2*M_PI/n, fi1=0;//M_PI/2;
  for (int i=0; i<n; i++)
    cnt->add(Vector_2(center[0]+d*cos(fi1+i*dfi), center[1]+d*sin(fi1+i*dfi)));
  return new Polygon_2(cnt);
}

int build_orth_basis(Vector_3 *in, int innum, Vector_3 *out){
  for(int i=0;i<innum;i++){
    if(in[i].norm()<=VEC_ZERO)
      return 0;
    in[i].normalize();
  }

  if(innum==0){
    for(int i=0;i<3;i++){
      out[i]=Vector_3();
      out[i][i]=1;
    }
  }
  if(innum==1){
    for(int i=0;i<3;i++){
      int i1=(i+1)%3;
      int i2=(i+2)%3;
      out[0][i] = in[0][i1]-in[0][i2];
    }
    out[0].normalize();
    out[1]=in[0]%out[0];
  }
  if(innum==2){
    out[0]=in[0]%in[1];
    if(out[0].norm()<=VEC_ZERO)
      return 0;
  }
  return 1;
}
