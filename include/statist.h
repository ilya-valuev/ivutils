/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.14 $
 *   $Date: 2015/11/09 18:16:03 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/statist.h,v 1.14 2015/11/09 18:16:03 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/statist.h,v $
$Revision: 1.14 $
$Author: valuev $
$Date: 2015/11/09 18:16:03 $
*/
/*e****************************************************************************
 * $Log: statist.h,v $
 * Revision 1.14  2015/11/09 18:16:03  valuev
 * added run average counter
 *
 * Revision 1.13  2015/11/05 11:22:28  valuev
 * fixed trajectory saving (regular sequence)
 *
 * Revision 1.12  2014/07/16 10:32:11  valuev
 * quantum sums for AWPMC
 *
 * Revision 1.11  2013/08/16 11:32:45  valuev
 * compiled tcpengine
 *
 * Revision 1.10  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.9  2009/06/10 20:53:33  valuev
 * updated splindes, trajReader
 *
 * Revision 1.8  2009/02/10 14:20:45  valuev
 * sync with FDTD project
 *
 * Revision 1.4  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.7  2008/03/18 17:17:21  valuev
 * corrected PBC in microfield calculations
 *
 * Revision 1.6  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.5  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.2  2006/08/08 13:02:04  valuev
 * Added geometry
 *
 * Revision 1.3  2006/03/14 10:32:17  valuev
 * Added SetControls support for many components,
 * improved tcpenfine, added GRASP interface
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
#ifndef __STATISTICS_H
#define __STATISTICS_H

#include <math.h>
#include "math_utils.h"  // accdiv
#include "common.h"

enum{
 VC_FFT,
 VC_DIRECT,
 VC_STRAIT
};



class Statistics{
protected:
  virtual double pnext(int dn=1){
   s+=dn*cur;
   s2+=dn*cur*cur;
   nsteps+=dn;
   return s/nsteps;
  }
public:
  double cur;
  double *value;
  double s; //e< sum
  double s2; //e< sum of squares
  long nsteps;

  virtual void clear(){
   s=s2=0.;
   nsteps=0;
  }

  Statistics(){
   cur=0.;
   value=&cur;
   clear();
  }

  virtual Statistics &operator=(const Statistics &a){
   value=a.value;
   s=a.s;
   s2=a.s2;
   nsteps=a.nsteps;
   return *this;
  }

  virtual Statistics &operator+=(const Statistics &a){
   s+=a.s;
   s2+=a.s2;
   nsteps+=a.nsteps;
   return *this;
  }

  virtual Statistics &operator/=(double r){
   s/=r;
   s2/=r*r;
   return *this;
  }
  virtual Statistics &operator*=(double r){
   s*=r;
   s2*=r*r;
   return *this;
  }
  Statistics(double *p){
   cur=0.;
   value=p;
   clear();
  }

  void attach(double *p=NULL){
    if(!p)value=&cur;
    else value=p;
  }


  double next(int dn=1){
   cur=value[0];
   return pnext(dn);
  }
  double next(double s, int dn=1){
   cur=s;
   return pnext(dn);
  }
  double prev(int dn=1){
   return pnext(dn);
  }
  double av(){
   if(nsteps)return s/nsteps;
   else return 0.;
  }
  double av2(){
   if(nsteps)return s2/nsteps;
   else return 0.;
  }
  double dev(){
   if(nsteps)return sqrt(fabs((s2-s*s/nsteps)/nsteps));
   else return 0.;
  }
  double rel(){
   if(nsteps)return dev()/fabs(s/nsteps);
   else return 0.;
  }
  double av(int n){
   return s/n;
  }
  double av2(int n){
   return s2/n;
  }
  double dev(int n){
   return sqrt(fabs((s2-s*s/n)/n));
  }
  double rel(int n){
   return dev(n)/fabs(s/n);
  }
};

//e this is correct statistics
class Statistics2{
protected:
  virtual double pnext(int dn=1){
    int i;
    for(i=0;i<dn;i++){
      if(!nsteps){
        s=cur;
        s2=cur*cur;
        sn=sxn=sn2=0;
        nsteps++;
      }
      else{
        double icur=nsteps;
        nsteps++;
        s+=(cur-s)/nsteps;
        s2+=(cur*cur-s2)/nsteps;
        sn+=(icur-sn)/nsteps;
        sxn+=(cur*icur-sxn)/nsteps;
        sn2+=(icur*icur-sn2)/nsteps;
      } 
    }
    return s;
  }
  virtual double ipnext(double icur){
    if(!nsteps){
      s=cur;
      s2=cur*cur;
      sn=icur;
      sxn=cur*icur;
      sn2=icur*icur;
      nsteps++;
    }
    else{
      nsteps++;
      s+=(cur-s)/nsteps;
      s2+=(cur*cur-s2)/nsteps;
      sn+=(icur-sn)/nsteps;
      sxn+=(cur*icur-sxn)/nsteps;
      sn2+=(icur*icur-sn2)/nsteps;
    }    
    return s;
  }
  
public:
  double cur;
  double *value;
  double s; //e< sum of x/n
  double s2; //e< sum of x squares/n
  double sxn; //e< sum of x*i /n
  double sn2; //e< sum of i*i /n
  double sn;  //e< sum of i/n
  long nsteps;

  virtual void clear(){
   sn=sxn=sn2=s=s2=0.;
   nsteps=0;
  }

  Statistics2(){
   cur=0.;
   value=&cur;
   clear();
  }

  virtual Statistics2 &operator=(const Statistics2 &a){
    cur=a.cur;
    value=a.value;
    s=a.s;
    s2=a.s2;
    sxn=a.sxn;
    sn2=a.sn2;
    sn=a.sn;
    nsteps=a.nsteps;
    return *this;
  }

  virtual Statistics2 &operator+=(const Statistics2 &a){
    double ntot=nsteps+a.nsteps;
    double w1=nsteps/ntot;
    double w2=a.nsteps/ntot;
    s=s*w1+a.s*w2;
    s2=s2*w1+a.s2*w2;
    sxn=sxn*w1+a.sxn*w2;
    sn2=sn2*w1+a.sn2*w2;
    sn=sn*w1+a.sn*w2;
    nsteps+=a.nsteps;
    return *this;
  }

  virtual Statistics2 &operator/=(double r){
   s/=r;
   sxn/=r;
   s2/=r*r;
   return *this;
  }
  virtual Statistics2 &operator*=(double r){
   s*=r;
   s2*=r*r;
   sxn*=r;
   return *this;
  }
  Statistics2(double *p){
   cur=0.;
   value=p;
   clear();
  }

  void attach(double *p=NULL){
    if(!p)value=&cur;
    else value=p;
  }


  double next(int dn=1){
   cur=value[0];
   return pnext(dn);
  }
  double next(double s, int dn=1){
   cur=s;
   return pnext(dn);
  }
  
  double nextxy(double x,double y){
   cur=y;
   return ipnext(x);
  }

  double prev(int dn=1){
   return pnext(dn);
  }
  double av(){
   return s;
  }
  double av2(){
   return s2;
  }

  double dev(){
    return sqrt(fabs(s2-s*s));
  }

  double rel(){
   if(nsteps)return dev()/fabs(s);
   else return 0.;
  }
  double av(int n){
   return s*nsteps/n;
  }
  double av2(int n){
   return s2*nsteps/n;
  }
  double dev(int n){
   return sqrt(fabs(s2-s*s)*nsteps/n);
  }
  double rel(int n){
   return dev(n)/fabs(s*nsteps/n);
  }
  //e returns linear trend x = a*i+b
  //e in a and b
  //e and deviation from it as return value
  double lintrend(double *a=NULL, double *b=NULL){
    if(nsteps<2)return -1.;
    //double sn=(nsteps-1)/2.;
    double devx2=s2-s*s;
    double devn2=sn2-sn*sn;
    double ta=0.;
    if(devx2>1e-32){
      ta=(sxn-sn*s)/devn2;
    }
    double tb=s-ta*sn;
    double tdev=devx2-ta*ta*devn2-2*ta*(sxn-s*sn);
    if(a)*a=ta;
    if(b)*b=tb;
    return sqrt(fabs(tdev));
  }

 

};

typedef Statistics2 *Statistics2P;


class RunAverage: public Statistics2{
protected:
  double pnext(int dn=1);


  double *array;
public:
  long cur_n;
  long window;
  
  RunAverage(long window_=1):Statistics2(){
    window=window_;
    cur_n=0;
    array= new double[window];
    if(!array)msg_error("RunAver: allocation error!\n");
  }

  RunAverage(double *p, long w=1):Statistics2(p){
    window=w;
    cur_n=0;
    array=new double[window];
    if(!array)msg_error("RunAver: allocation error!\n");
  }


  virtual ~RunAverage(){
    delete[] array;
  }

  int set_window(long w);

  int pushout(int dn=1);

  void clear(){
    cur_n=0;
    Statistics2::clear();
  }
  //e average
  double rav() const;
  //e average square
  double rav2() const;
   //e running minimum
  double rmin() const;
  //e running maximum
  double rmax() const;

  //e standart deviation
  double rdev() const{
    if(cur_n){
      double ra=rav();
      double ra2=rav2();
      return sqrt(fabs( (ra2 -   ra*ra)));
    }
    else return 0.;
  }

  //e deviation relative to average
  double rrel() const {
    if(cur_n)return rdev()/fabs(rav());
    else return 0.;
  }



  Statistics2 &operator=(Statistics2 a){
    eprintf("RunAverage: operation '=' not allowed !\n");
    return *this;
  }

  Statistics2 &operator+=(Statistics2 a){
    eprintf("RunAverage: operation '+=' not allowed !\n");
    return *this;
  }

  Statistics2 &operator/=(double r);

  Statistics2 &operator*=(double r);


};

typedef RunAverage *RunAverageP;


void Smooth(TableFunction *func,int window);




class Distribution{
public:
 TableFunction distr;
 Statistics2 prop;
 double value;
 realtype *arr;

public:
 realtype x1,x2;
 realtype dx;
 realtype norm;
 realtype norm_r2;
 realtype norma;
 realtype norm_r2a;

 int n;
 unsigned long count;
 unsigned long allcount;
 Distribution();
 Distribution(realtype x1,realtype x2, int n); //:prop(&value);
 
 //e can be initialized only once
 //e returns >0 if OK
 //e <0 if error
 int init(realtype x1,realtype x2, int n);
 
 ~Distribution(){ 
   if(arr)delete[] arr; 
 }
 
 realtype operator()(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   return distr(x);
 }
 realtype nrm(realtype x){
   if(!count)return 0;
   return distr(x)/norm;
 }
 realtype nrm_r2(realtype x){
   if(!count)return 0;
   return distr(x)*x*x/norm_r2;
 }

 realtype nrma(realtype x){
   if(!count)return 0;
   return distr(x)/norma;
 }
 realtype nrma_r2(realtype x){
   if(!count)return 0;
   return distr(x)*x*x/norm_r2a;
 }
 ///\en Gets r2 norm calculated up to limiting argument value
 realtype get_norm_r2(realtype xlim, realtype dx_=-1.) {
   if (!count)return 0.;
   double mult = 0, x;
   if (dx<0)
     dx = (x2 - x1) / (3 * (n - 1));
   for (x = x1; x<xlim; x += dx)
     mult += distr(x)*x*x*dx;
   return isnan(1./mult) ?  1. : mult ;
 }


 int point(realtype x,realtype  val);

 ///\en same as above but with 3d weight 1/(4*pi*x*x*dx)  
 int point_r2(realtype x, realtype  val) {
   double k = 4 * M_PI*dx;
   double f = fabs(x) < dx ? 1. / (k * dx*dx / 3.) : 1. / (k*x*x);
   return point(x, f);
 }

 realtype normalize(realtype norm=1.);
 realtype normalize_r2(realtype norm=1.,realtype dx=-1.);

 //e writes distribution to file
 //e when normalized =0 does not normalize
 //e                 =1 normalize to thenorm
 //e                 =2 normalize to thenorm with points outside the interval
 //e                 =3 normalize to thenorm integrating with x*x
 //e                 =4 normalize to thenorm integrating with x*x with points outside the interval
 void write(char *file,int normalized=0, realtype thenorm=1.){
   realtype scale=thenorm;
   if(normalized==1)scale/=norm;
   else if(normalized==2)scale/=norma;
   else if(normalized==3)scale/=norm_r2;
   else if(normalized==4)scale/=norm_r2a;
   distr.write(file,scale);
 }
 realtype av(){
  return (realtype)prop.av();
 }
 realtype av2(){
  return (realtype)prop.av2();
 }
 realtype dev(){
  return (realtype)prop.dev();
 }
 void clear();

 //void *set(realtype (*func)(realtype x));

 TableFunction *get_func(){
   return &distr;
 }
};


enum {CORR_COMPLEX=1, CORR_REAL=0};
enum {CORR_DIRECT=0x10,CORR_FFT=0x20};

class Correlation{
public:

 realtype *ini; // initial data
 realtype *inii;
 TableFunction fini;
 TableFunction finii; // im part


 TableFunction corr; // correlation
 TableFunction corri; // im part
 TableFunction four_sq;
 TableFunction four_sqi;
 realtype *arr;
 realtype *arri;

 int n;
 int status;
 realtype maxtime;
 realtype sc;

public:
 int type;
 int sub_mean;
 //int method;

 Correlation();
 Correlation(int n, realtype maxtime=0., int ctype=CORR_REAL);
 void init(int n, realtype maxtime, int ctype=CORR_REAL);
 ~Correlation(){
  if(arr)delete[] arr;
  if(arri)delete[] arri;
  if(ini)delete[] ini;
  if(inii)delete[] inii;
 }
 realtype operator()(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   return corr(x);
 }
 realtype r(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   return corr(x);
 }
 realtype i(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   if(type==CORR_COMPLEX)return corri(x);
   else return 0.;
 }

 TableFunction& func_r(){
   return corr;
 }

 TableFunction& func_i(){
   return corri;
 }

 TableFunction& func_fsqr(){
   return four_sq;
 }

 TableFunction& func_fsqi(){
   return four_sqi;
 }


 realtype data(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   return fini(x);
 }
 realtype data_r(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   return fini(x);
 }
 realtype data_i(realtype x){
   //printf("pointer %ld  %ld\n",(long)distr.yy,(long)arr);
   if(type==CORR_COMPLEX)return finii(x);
   else return 0.;
 }


 void insert_begin(int steps, realtype factor=0.){
   if(factor==0.){
     sc=(realtype)sqrt((realtype)n/steps);
   }
   else sc=factor;

   status=0;
   //corr.insert_begin(steps);
   fini.insert_begin(steps);
   if(type==CORR_COMPLEX){
    //corri.insert_begin(steps);
    finii.insert_begin(steps);
   }
 }

 realtype insert_next(realtype x,realtype y=0.){

   if(type==CORR_COMPLEX){
    //corri.insert_next(y);
    finii.insert_next(y*sc);
   }
   return fini.insert_next(x*sc);
   //corr.insert_next(x);
 }

 void insert_end();

 realtype  *calculate();
 void direct();
 void strait();
 realtype normalize();
 void set_mean(realtype meanv=0., realtype meanvi=0.);
 void set_copy();
};

realtype *Calculate_forth(int n, Correlation *corrset, int num);
realtype *Calculate_back(realtype *arr, int n);


# define DIFFSZ  // allow different sizes for c1 and c2 in cross correlations

void CrossCorrFFT(Correlation *res12, Correlation *res21, Correlation *c1, Correlation *c2);
void CrossCorrDirect(Correlation *res12, Correlation *res21, Correlation *c1, Correlation *c2);
void CrossCorrStrait(Correlation *res12, Correlation *res21, Correlation *c1, Correlation *c2);

//e Sequence linear interpolator
class SeqInterp{
  double t0, dt;
  double tw;
public:
  //e How to use it:
  //e 1. constuct SeqInterp and insert A(t0)
  //e 2. store A0=A(t0)
  SeqInterp(double st0=0, double sdt=1.){
    init(st0,sdt);
  }
  int init(double st0, double sdt){
    t0=st0;
    dt=sdt;
    tw=t0+dt;
    return 1;
  }
  //e 3. for each A(ti)
  //e call n=next(ti,&cl1,&dcl1)
  //e then for i=0;i<n;i++
  //e c=cl+dcl*i
  //e insert c*A0+(1-c)*A(ti)
  //e store A0=A(ti)
  int next(double ti, double *cl, double *dcl){
    if(ti<tw){
      t0=ti;
      return 0;
    }
    int nu=(int)(accdiv<double>(ti-tw,dt))+1;
    int nl=(int)(accdiv<double>(tw-t0,dt));
    if(acccomp(tw-nl*dt,t0))// exact ?
      nl--;
    double l=ti-t0;
    *cl=1.-(tw-nl*dt-t0)/dt;
    *dcl=-dt/l;
    t0=ti;
    tw+=nu*dt;
    return nu+nl;
    /*double l=ti-t0, d1=tw-t0;
    *cl=1.-d1/l;
    *dcl=-dt/l;
    int n=(int)((ti-tw)/dt);
    tw+=(n+1)*dt;
    t0=ti;
    return n+1;*/
  }
};


# endif













