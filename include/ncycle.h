#ifndef _NCYCLE_H
#define _NCYCLE_H

# include "common.h"

class CycleParam{
  double p1;
  double p2;
  double dp;
  int i;
  int n;
  CycleParam *nextp;
  char name[250];
  char digfmt[10];
public:
  int init(char *spec,double start, double end, double step);
  int init(char *spec,char *range, char *delim=",");
  
  CycleParam(char *spec,double start, double end, double step){
    init(spec,start,end,step);
  }  

  CycleParam(char *spec,double value){
    init(spec,value,value,1.);
  }
  CycleParam(char *spec,char *range, char *delim=","){
    init(spec,range,delim);
  }

  CycleParam(){
    n=0;
    i=0;
    p1=p2=dp=0;
  }

  double value(){
    if(n>1)return p1+(p2-p1)*i/(n-1);
    else return p1;
  }  
  char *format(){
    return digfmt;
  }
  char *getname(){
    return name;
  }

  int nsteps(){
    return n;
  }

  friend class Cycle;
};

typedef CycleParam *CycleParamP;


class Cycle{
  CycleParam** p;
  int maxnp;
  int np;

  int init(int n);
  char buffer[250];

public:
  Cycle(int n){
    if(!init(n))fatal_error("MAE at Cycle allocation!\n");
  }

  int Add(CycleParam *par){
    p[np]=par;
    np++;
    return np;
  }

  int Next();
  int Current();
  int Total(int lev=-1);
  virtual char *SpecString();
 
  // does not delete parameters !!!
  virtual ~Cycle(){
    delete p;
  }

};





#endif
