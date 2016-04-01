# ifndef SECTION_H
# define SECTION_H

# include "curve2d.h"

# define INI_NPOINTS 10


class SectionCurve: public BoxedCurve{
protected:
  Pool pval;
  int np;
  float lastval;
  void to_integral(int num);

  TableFunction fvl;

  int end_input(){
    if(!input)return 0;
    fvl=TableFunction(np,ampl);
    BoxedCurve::end_input();
    return 1;
  }

public:

  float sum;
  float *ampl;

  float Integral();

  SectionCurve():BoxedCurve(),pval((void **)&ampl,(int)sizeof(float)){
    //printf("scccon ");
  }

  int init(float x1,float y1, float x2, float y2, float bw){
    pval.clear();
    pval.size1=INI_NPOINTS;
    pval.size_add=D_NPOINTS;
    np=0;
    lastval=0.;
    sum=0;
    return BoxedCurve::init(INI_NPOINTS,x1,y1,x2,y2,bw);
  }

  SectionCurve(float x1,float y1, float x2, float y2, float bw):
    pval((void **)&ampl,(int)sizeof(float)){
      init(x1,y1,x2,y2,bw);
  }


  ~SectionCurve(){
    if(ampl)free(ampl);
  }

  int AddPointVal(float x, float y, float val);
  int AddPointRequest(float x, float y);
  int FillRequested(float val, int ind=-1); 

  void BeginRequests(){
    np=0;
    lastval=0.;
    sum=0.;
  }

  int AddPoint(float x, float y){
    return AddPointVal(x,y,0.);
  }

  float value(int ind,float t=0.){
    if(ind<0 || ind>=np)return 0.;
    if(ind==np-1)return ampl[ind];
    return ampl[ind]*(1.-t)+ampl[ind+1]*t;
  }

  float value_l(float l);
  float value_lnorm(float t);

  int out(char *file);

};




# endif












