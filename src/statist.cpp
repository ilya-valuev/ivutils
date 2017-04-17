/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.6 $
 *   $Date: 2015/11/09 18:16:03 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/statist.cpp,v 1.6 2015/11/09 18:16:03 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/statist.cpp,v $
$Revision: 1.6 $
$Author: valuev $
$Date: 2015/11/09 18:16:03 $
*/
/*e****************************************************************************
 * $Log: statist.cpp,v $
 * Revision 1.6  2015/11/09 18:16:03  valuev
 * added run average counter
 *
 * Revision 1.5  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.4  2007/06/15 13:10:10  valuev
 * modified geometry for oblique incidence
 *
 * Revision 1.3  2006/10/27 20:41:02  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.3  2006/10/09 12:58:56  bogomolov
 * patch for compilling under unix
 *
 * Revision 1.2  2006/03/14 10:32:17  valuev
 * Added SetControls support for many components,
 * improved tcpenfine, added GRASP interface
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:16  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/
# include "statist.h"             
# include "four.h" 


# if 0   //  ndef FUCK
double log2(double x);


double log2(double x){
  return log(x)/log(2);
}

# endif

//double Statistics::cur=0.;

// centered smoothing
void Smooth(TableFunction *func,int window){
  RunAverage s(NULL,window);
  int i;
  int right=window/2, top=func->n+right;

  for(i=0;i<right;i++)s.next(func->y(i));  // right half-window
  for(;i<func->n;i++){
    s.next(func->y(i));
    func->y(i-right)=s.rav();
  }
  for(;i<top;i++){
    s.pushout();
    func->y(i-right)=s.rav();
  }
}


Distribution::Distribution(realtype x1,realtype x2, int np):prop(&value),arr(NULL){
  init(x1,x2,np);
}

Distribution::Distribution():prop(&value){
  arr=NULL;
  n=0;
}

int Distribution::init(realtype x1,realtype x2, int np){
  if(arr){
    delete [] arr;
    arr=NULL;
  }
  arr=new realtype[np];
  n=np;
  if(!arr){
    msg_error("Distribution: can not allocate storage array!\n");
    return -2;
  }
  int i;
  for(i=0;i<np;i++)arr[i]=0;
  TableFunction t_distr(n,arr);
  distr=t_distr;
  //distr=TableFunction(n,arr);
  distr.xscale(x1,x2);
  // clear();
  //printf("constr --- %f  %f\n",distr(x1),distr(0));
  this->x1=x1;
  this->x2=x2;
  count=0;
  allcount=0;
  norm=norm_r2=0;
  norma=norm_r2a=0;
  return 1;
}

void Distribution::clear(){
  if(!arr)return;
  int i;
  for(i=0;i<n;i++)arr[i]=0;
  count=0;
  allcount=0;
  norm=norm_r2=0;
  prop.clear();
}

int Distribution::point(realtype x,realtype  val){
  allcount++;
  
  realtype dx=(x2-x1)/(n-1);
  realtype x0=x1-0.5*dx;
  int ind=(int)( (x-x0)/dx);
  if(ind==0 || ind==n-1)val*=2; // half-interval correction
  double d1=dx*val, d2=d1*x*x;
  norma+=d1;
  norm_r2a+=d2;
  if(x<x1 || x>x2){
    return -1; //prop.av();
  }
  arr[ind]+=val;
  count++;
  norm+=d1;
  norm_r2+=d2;
  value=x;
  prop.next();
  return ind;
}

realtype Distribution::normalize(realtype new_norm){
  if(!count)return 0.;
  int i;
  double mult=new_norm/norm;
  for(i=0;i<n;i++)arr[i]*=mult;
  norm=new_norm;
  norm_r2*=mult;
  //prop.s*=mult;
  //prop.s2*=(mult*mult);
  return prop.av();
}



realtype Distribution::normalize_r2(realtype new_norm_r2, realtype dx){
  if(!count)return 0.;
  int i;
  double mult=0,x;
  if(dx<0)dx=(x2-x1)/(3*(n-1));
  for(x=x1;x<x2;x+=dx)mult+=distr(x)*x*x*dx;

  mult=new_norm_r2/mult;
  for(i=0;i<n;i++)arr[i]*=mult;
  norm*=mult;
  norm_r2*=mult;
  //prop.s*=mult;
  //prop.s2*=(mult*mult);
  return prop.av();
}


# ifdef UNIX
Correlation::Correlation(int n, realtype smaxtime,int ctype){
# else
Correlation::Correlation(int n, realtype smaxtime,int ctype){
# endif
  ini=NULL;
  inii=NULL;
  arr=NULL;
  arri=NULL;
  if(smaxtime<=0.)smaxtime=(realtype)n;
  init(n,smaxtime,ctype);
}


Correlation::Correlation(){
  ini=NULL;
  inii=NULL;
  arr=NULL;
  arri=NULL;
  n=0;
  type=CORR_REAL;
  //method=CORR_FFT;
  sub_mean=0;
}

# ifndef FU
int rint(double x){
  int fl=(int)floor(x);
  if(x-(double)fl > 0.5)return fl+1;
  else return fl;
}
# endif


void Correlation::init(int np, realtype smaxtime,int ctype){
  type=ctype;


  if(this->arr)delete[] this->arr;
  if(this->arri)delete[] this->arri;
  if(this->ini)delete[] this->ini;
  if(this->inii)delete[] this->inii;
  

  this->arr=new realtype[2*np];
  n=np;
  if(!this->arr)serror("Correlation: can not allocate storage array!\n");
  //if(type==CORR_COMPLEX){
  this->arri=new realtype[2*np];
  if(!this->arri)serror("Correlation: can not allocate storage array!\n");
  //}



  ini=new realtype[np];
  if(!ini)serror("Correlation: can not allocate data array!\n");

  if(type==CORR_COMPLEX){
   inii=new realtype[np];
   if(!inii)serror("Correlation: can not allocate data array!\n");
  }

  int i;
  for(i=0;i<2*np;i++)arr[i]=0;
  TableFunction t_corr(n,arr),t_four_sq(n,arr+n),t_four_sqi(n,arri+n),t_fini(n,ini);
  corr=t_corr;
  four_sq=t_four_sq;
  four_sqi=t_four_sqi;
  fini=t_fini;
  //corr=TableFunction::TableFunction(n,arr);
  //four_sq=TableFunction::TableFunction(n,arr+n);
  //four_sqi=TableFunction::TableFunction(n,arri+n);
  //fini=TableFunction::TableFunction(n,ini);

  maxtime = smaxtime;
  corr.xscale(0.,maxtime);
  fini.xscale(0.,maxtime);
  four_sq.xscale(0.,2.*M_PI*n/maxtime);
  four_sqi.xscale(0.,2.*M_PI*n/maxtime);

  TableFunction t_corri(n,arri);
  corri=t_corri;
  //corri=TableFunction(n,arri);
  corri.xscale(0.,maxtime);
  if(type==CORR_COMPLEX){
   TableFunction t_finii(n,inii);
   finii=t_finii;
   //finii=TableFunction(n,inii);
   finii.xscale(0.,maxtime);
  }

  status=0;
  sub_mean=0;
}


realtype *Calculate_forth(int n,Correlation *corrset, int num){
  int i,l;

  realtype *cur_arr=(realtype *)malloc(4*n*sizeof(realtype));
  //new realtype[4*n];
  realtype *tmp_arr= new realtype[8*n];
  if(!cur_arr || !tmp_arr)serror("Correlation.calculate: memory allocation error.\n");


  for(i=0;i<4*n;i++)tmp_arr[i]=0;

  double sum=0.;
  for(l=0;l<num;l++){

    if(corrset[l].sub_mean)corrset[l].set_mean();
    else corrset[l].set_copy();


    for(i=0;i<n;i++){
      cur_arr[2*i]=corrset[l].arr[i];
      cur_arr[2*(2*n-i-1)]= 0.; // cur_arr[2*i];
      if(corrset[l].type==CORR_COMPLEX){
       cur_arr[2*i+1]=corrset[l].arri[i];
       cur_arr[2*(2*n-i-1)+1]=0.; //cur_arr[2*i+1];
      }
      else{
       cur_arr[2*(2*n-i-1)+1]=0.;
       cur_arr[2*i+1]=0.;
      }
    }

    //four(cur_arr,2*n,1);
    four(cur_arr,n,1);

    for(i=0;i<n;i++){
      cur_arr[2*i]/=n; // scaling after fourier transform
      cur_arr[2*i+1]/=n;
      corrset[l].arr[i]=cur_arr[2*i];  // copying fourier transform
      corrset[l].arri[i]=cur_arr[2*i+1];
      // summing fourier transform square
      corrset[l].arr[n+i]=cur_arr[2*i]*cur_arr[2*i]+cur_arr[2*i+1]*cur_arr[2*i+1];
      corrset[l].arri[n+i]=0.;
      tmp_arr[2*i]+=corrset[l].arr[n+i];
      sum+=corrset[l].arr[n+i];
    }


  }
  sum/=(n*num);
  for(i=0;i<n;i++){
    tmp_arr[2*i]/=num;
    if(corrset[0].sub_mean)
      tmp_arr[2*i]-=sum;
    tmp_arr[2*i+1]=0.;
  }
  for(i=0;i<n;i++){
    tmp_arr[2*n+i]=tmp_arr[2*i];//copying fourier transform square

  }

  free(cur_arr);
  return tmp_arr;
}

realtype *Calculate_back(realtype *fft_sq, int n){
   //four(fft_sq,2*n,-1);
   //return fft_sq+4*n;
   /*int i;
   realtype sc=1./(n*n);
   for(i=0;i<2*n;i++){
    fft_sq[i]*=sc;
   } */
   four(fft_sq,n,-1);
   return fft_sq+2*n;
}


void Correlation::set_mean(realtype meanv, realtype meanvi){
 realtype sum=0.;
 int i;

 for(i=0;i<n;i++)sum+=ini[i];
 sum/=n;
 for(i=0;i<n;i++)arr[i]=ini[i]+meanv-sum; // substracting average

 if(type==CORR_COMPLEX){
  sum=0.;
  for(i=0;i<n;i++)sum+=inii[i];
  sum/=n;
  for(i=0;i<n;i++)arri[i]=inii[i]+meanvi-sum; // substracting average
 }
}

void Correlation::set_copy(){
 int i;
 for(i=0;i<n;i++)arr[i]=ini[i];
 if(type==CORR_COMPLEX){
  for(i=0;i<n;i++)arri[i]=inii[i];
 }
}


# ifndef DIFFSZ

void CrossCorrDirect(Correlation *res12, Correlation *res21, Correlation *c1, Correlation *c2){
  if(res12->n!=c1->n || res21->n!=c1->n || res12->n!=c2->n || c1->n!=c2->n){
   serror("CrossCorrDirect: sizes do not match!\n");
  }
  int n=res12->n;
  realtype sum12, sumi12, sum21, sumi21 ;
  int i,j;

  if(res12->sub_mean){
   c1->set_mean();
   c2->set_mean();
  }
  else{
   c1->set_copy();
   c2->set_copy();
  }


  for(i=0;i<n;i++){
   sum12=0.;
   sumi12=0.;
   sum21=0.;
   sumi21=0.;
   for(j=0;j+i<n;j++){      // c1*(j)mal c2(j+i) , c2*(j) mal c1(j+i)

    sum12+=c1->arr[j]*c2->arr[j+i];
    sum21+=c2->arr[j]*c1->arr[j+i];
    if(res12->type==CORR_COMPLEX){

     sum12+=c1->arri[j]*c2->arri[j+i];
     sum21+=c2->arri[j]*c1->arri[j+i];


     sumi12-=c1->arri[j]*c2->arr[j+i]-c1->arr[j]*c2->arri[j+i];
     sumi21-=c2->arri[j]*c1->arr[j+i]-c2->arr[j]*c1->arri[j+i];
    }
   }
   res12->arr[i]=sum12/(n-i);
   res21->arr[i]=sum21/(n-i);

   if(res12->type==CORR_COMPLEX){
    res12->arri[i]=sumi12/(n-i);
    res21->arri[i]=sumi21/(n-i);
   }
  }
  res12->status=1;
  res21->status=1;
  c1->direct();
  c2->direct();
}

# else



void CrossCorrDirect(Correlation *res12, Correlation *res21, Correlation *c1, Correlation *c2){
  if(res12->n!=res21->n){
   serror("CrossCorrDirect: sizes do not match!\n");
  }
  if(fabs(res12->maxtime-res21->maxtime)> 1e-10){
   serror("CrossCorrDirect: time scales do not match!\n");
  }
  int n=res12->n;
  realtype sum12, sumi12, sum21, sumi21 ;
  int i,j;

  if(res12->sub_mean){
   c1->set_mean();
   c2->set_mean();
  }
  else{
   c1->set_copy();
   c2->set_copy();
  }

  realtype tj,tij;
  realtype tscale=res12->maxtime;


  for(i=0;i<n;i++){
   sum12=0.;
   sumi12=0.;
   sum21=0.;
   sumi21=0.;
   for(j=0;j+i<n;j++){      // c1*(j)mal c2(j+i) , c2*(j) mal c1(j+i)

    tj=j*tscale/(n-1);
    tij=(j+i)*tscale/(n-1);

    sum12+=c1->r(tj)*c2->r(tij);
    sum21+=c2->r(tj)*c1->r(tij);
    if(res12->type==CORR_COMPLEX){

     sum12+=c1->i(tj)*c2->i(tij);
     sum21+=c2->i(tj)*c1->i(tij);


     sumi12-= c1->i(tj)*c2->r(tij) - c1->r(tj)*c2->i(tij);
     sumi21-= c2->i(tj)*c1->r(tij) - c2->r(tj)*c1->i(tij);
    }
   }
   res12->arr[i]=sum12/(n-i);
   res21->arr[i]=sum21/(n-i);

   if(res12->type==CORR_COMPLEX){
    res12->arri[i]=sumi12/(n-i);
    res21->arri[i]=sumi21/(n-i);
   }
  }
  res12->status=1;
  res21->status=1;
  c1->direct();
  c2->direct();
}

# endif

# ifndef DIFFSZ


void CrossCorrStrait(Correlation *res12, Correlation *res21,
		     Correlation *c1, Correlation *c2){

 if(res12->n!=c1->n || res21->n!=c1->n || res12->n!=c2->n || c1->n!=c2->n){
  serror("CrossCorrStrait: sizes do not match!\n");
 }
 int n=res12->n;
 int i;

 if(res12->sub_mean){
   c1->set_mean();
   c2->set_mean();
 }
 else{
   c1->set_copy();
   c2->set_copy();
 }

 if(res12->type==CORR_COMPLEX){
  realtype x1=c1->arr[0];
  realtype y1=c1->arri[0];
  realtype x2=c2->arr[0];
  realtype y2=c2->arri[0];
  for(i=0;i<n;i++){
   res12->arr[i] =(c2->arri[i]*x1+c2->arri[i]*y1);  // c1*(0)mal c2(t)
   res12->arri[i]=(c2->arri[i]*x1-c2->arr[i]*y1);

   res21->arr[i] =c1->arr[i]*x2+c1->arri[i]*y2;  // c2*(0)mal c1(t)
   res21->arri[i]=c1->arri[i]*x2-c1->arr[i]*y2;
  }
 }
 else{
  realtype x1=c1->arr[0];
  realtype x2=c2->arr[0];
  for(i=0;i<n;i++){
   res12->arr[i]=c2->arr[i]*x1;
   res21->arr[i]=c1->arr[i]*x2;
  }
 }


 /*
 if(res12->type==CORR_COMPLEX){
  realtype x1=c1->ini[0];
  realtype y1=c1->inii[0];
  realtype x2=c2->ini[0];
  realtype y2=c2->inii[0];
  for(i=0;i<n;i++){
   res12->arr[i] =c2->ini[i]*x1+c2->inii[i]*y1;  // c1*(0)mal c2(t)
   res12->arri[i]=c2->inii[i]*x1-c2->ini[i]*y1;

   res21->arr[i] =c1->ini[i]*x2+c1->inii[i]*y2;  // c2*(0)mal c1(t)
   res21->arri[i]=c1->inii[i]*x2-c1->ini[i]*y2;
  }
 }
 else{
  realtype x1=c1->ini[0];
  realtype x2=c2->ini[0];
  for(i=0;i<n;i++){
   res12->arr[i]=c2->ini[i]*x1;
   res21->arr[i]=c1->ini[i]*x2;
  }
 } */
 c1->strait();
 c2->strait();
}

# else


void CrossCorrStrait(Correlation *res12, Correlation *res21,
		     Correlation *c1, Correlation *c2){

 if(res12->n!=res21->n){
   serror("CrossCorrStrait: sizes do not match!\n");
 }
 if(fabs(res12->maxtime-res21->maxtime)> 1e-10){
   serror("CrossCorrStrait: time scales do not match!\n");
 }
 realtype ti;
 realtype tscale=res12->maxtime;

 int n=res12->n;
 int i;

 if(res12->sub_mean){
   c1->set_mean();
   c2->set_mean();
 }
 else{
   c1->set_copy();
   c2->set_copy();
 }

 if(res12->type==CORR_COMPLEX){
  realtype x1=c1->r(0.);
  realtype y1=c1->i(0.);
  realtype x2=c2->r(0.);
  realtype y2=c2->i(0.);
  for(i=0;i<n;i++){
   ti=i*tscale/(n-1);

   res12->arr[i] =c2->r(ti)*x1 + c2->i(ti)*y1;  // c1*(0)mal c2(t)
   res12->arri[i]=c2->i(ti)*x1 - c2->r(ti)*y1;

   res21->arr[i] =c1->r(ti)*x2 + c1->i(ti)*y2;  // c2*(0)mal c1(t)
   res21->arri[i]=c1->i(ti)*x2 - c1->r(ti)*y2;
  }
 }
 else{
  realtype x1=c1->r(0.);
  realtype x2=c2->r(0.);
  for(i=0;i<n;i++){
   ti=i*tscale/(n-1);
   res12->arr[i]=c2->r(ti)*x1;
   res21->arr[i]=c1->r(ti)*x2;
  }
 }


 /*
 if(res12->type==CORR_COMPLEX){
  realtype x1=c1->data_r(0.);
  realtype y1=c1->data_i(0.);
  realtype x2=c2->data_r(0.);
  realtype y2=c2->data_i(0.);
  for(i=0;i<n;i++){
   ti=i*tscale/(n-1);

   res12->arr[i] =c2->data_r(ti)*x1 + c2->data_i(ti)*y1;  // c1*(0)mal c2(t)
   res12->arri[i]=c2->data_i(ti)*x1 - c2->data_r(ti)*y1;

   res21->arr[i] =c1->data_r(ti)*x2 + c1->data_i(ti)*y2;  // c2*(0)mal c1(t)
   res21->arri[i]=c1->data_i(ti)*x2 - c1->data_r(ti)*y2;
  }
 }
 else{
  realtype x1=c1->data_r(0.);
  realtype x2=c2->data_r(0.);
  for(i=0;i<n;i++){
   ti=i*tscale/(n-1);
   res12->arr[i]=c2->data_r(ti)*x1;
   res21->arr[i]=c1->data_r(ti)*x2;
  }
 } */
 c1->strait();
 c2->strait();
}

# endif


# ifndef DIFFSZ


void CrossCorrFFT(Correlation *res12,Correlation *res21, Correlation *c1, Correlation *c2){
  if(res12->n!=c1->n || res21->n!=c1->n || res12->n!=c2->n || c1->n!=c2->n){
   serror("CrossCorrFFT: sizes do not match!\n");
  }
  int n=res12->n;
  if(fabs(log2(n)-rint(log2(n)))>=1e-5){
   printf("Correlation: number of points must be 2^n for FFT!\n"
          "(%d specified)\n",n );
   printf("Trying direct...\n");
   CrossCorrDirect(res12,res21,c1,c2);
  }

  realtype *tmp1=Calculate_forth(n,c1, 1);
  realtype *tmp2=Calculate_forth(n,c2, 1);

  realtype *tmpres12=new realtype[2*n];
  realtype *tmpres21=new realtype[2*n];
  if(!tmpres12 || !tmpres21)serror("CrossCorrFFT: memory allocation error.\n");

  int i;
  for(i=0;i<n;i++){
   res12->arr[n+i]=tmpres12[2*i]=c1->arr[i]*c2->arr[i]+c1->arri[i]*c2->arri[i];
   res12->arri[n+i]=tmpres12[2*i+1]=c1->arr[i]*c2->arri[i]-c1->arri[i]*c2->arr[i];
   res21->arr[n+i]=tmpres21[2*i]=tmpres12[2*i];
   res21->arri[n+i]=tmpres21[2*i+1]=-tmpres12[2*i+1]; // conjugate
  }

  Calculate_back(tmp1,n);
  Calculate_back(tmp2,n);
  Calculate_back(tmpres12,n);
  Calculate_back(tmpres21,n);

  for(i=0;i<n;i++){
   c1->arr[i]=tmp1[2*i];
   c1->arri[i]=tmp1[2*i+1];


   c2->arr[i]=tmp2[2*i];
   c2->arri[i]=tmp2[2*i+1];


   res12->arr[i]=tmpres12[2*i];
   res12->arri[i]=tmpres12[2*i+1];


   res21->arr[i]=tmpres21[2*i];
   res21->arri[i]=tmpres21[2*i+1];


  }
  c1->status=1;
  c2->status=1;
  res12->status=1;
  res21->status=1;

  delete [] tmp1;
  delete [] tmp2;
  delete [] tmpres12;
  delete [] tmpres21;
}

# else

void CrossCorrFFT(Correlation *res12,Correlation *res21, 
		  Correlation *c1, Correlation *c2){
  if(res12->n!=res21->n){
    serror("CrossCorrFFT: sizes do not match!\n");
  }
  if(fabs(c1->maxtime-c2->maxtime)> 1e-10){
    serror("CrossCorrFFT: time scales do not match!\n");
  }
  realtype ti;
  realtype tscale=c1->maxtime;
  

  int n=res12->n;
  int n1=c1->n;
  int n2=c2->n;

  if(fabs(log2(n) -rint(log2(n))) >=1e-5 ||
     fabs(log2(n1)-rint(log2(n1)))>=1e-5 ||
     fabs(log2(n2)-rint(log2(n2)))>=1e-5 ){
   printf("Correlation: number of points must be 2^n for FFT!\n"
          "(%d specified)\n",n );
   printf("Trying direct...\n");
   CrossCorrDirect(res12,res21,c1,c2);
   return;
  }


  realtype *tmp1=Calculate_forth(n1,c1, 1);
  realtype *tmp2=Calculate_forth(n2,c2, 1);

  realtype *tmpres12=new realtype[2*n];
  realtype *tmpres21=new realtype[2*n];
  if(!tmpres12 || !tmpres21)serror("CrossCorrFFT: memory allocation error.\n");

  int i;
  for(i=0;i<n;i++){  // multiplication
   ti=i*tscale/(n-1); 
   res12->arr[n+i]=tmpres12[2*i]=    c1->r(ti)*c2->r(ti) + c1->i(ti)*c2->i(ti);
   res12->arri[n+i]=tmpres12[2*i+1]= c1->r(ti)*c2->i(ti) - c1->i(ti)*c2->r(ti);
   res21->arr[n+i]=tmpres21[2*i]=tmpres12[2*i];
   res21->arri[n+i]=tmpres21[2*i+1]=-tmpres12[2*i+1]; // conjugate
  }

  Calculate_back(tmp1,n1);
  Calculate_back(tmp2,n2);
  Calculate_back(tmpres12,n);
  Calculate_back(tmpres21,n);

  for(i=0;i<n1;i++){
   c1->arr[i]=tmp1[2*i];
   c1->arri[i]=tmp1[2*i+1];
  }
  c1->status=1;

  for(i=0;i<n2;i++){
   c2->arr[i]=tmp2[2*i];
   c2->arri[i]=tmp2[2*i+1];
  }
  c2->status=1;

  for(i=0;i<n;i++){
   res12->arr[i]=tmpres12[2*i];
   res12->arri[i]=tmpres12[2*i+1];
   res12->status=1;

   res21->arr[i]=tmpres21[2*i];
   res21->arri[i]=tmpres21[2*i+1];
   res21->status=1;

  }

  delete [] tmp1;
  delete [] tmp2;
  delete [] tmpres12;
  delete [] tmpres21;
}

# endif


realtype *Correlation::calculate(){
  if(status==1){
    printf("Correlation: already performed.\n");
    return this->arr+n;
  }
  if(fabs(log2(n)-rint(log2(n)))>=1e-5){
   printf("Correlation: number of points must be 2^n for FFT!\n"
          "(%d specified)\n",n );
   printf("Trying direct...\n");
   direct();
   return arr+n;
  }

  realtype *tmp=Calculate_forth(this->n,this, 1);
  Calculate_back(tmp,n);

  int i;
  for(i=0;i<n;i++)arr[i]=tmp[2*i];
  //if(type==CORR_COMPLEX){
  for(i=0;i<n;i++)arri[i]=tmp[2*i+1];
  //}
  //for(i=n;i<2*n;i++)arr[i]=tmp[2*n+2*i];
  delete [] tmp;
  status=1;
  return arr+n;
}


void Correlation::direct(){
  realtype *tmp_arr, *tmp_arri=NULL;

  tmp_arr=new realtype[n];
  if(!tmp_arr)serror("Correlation.direct: can't allocate array!\n");

  if(type==CORR_COMPLEX){
   tmp_arri=new realtype[n];
   if(!tmp_arri)serror("Correlation.direct: can't allocate array!\n");
  }

  int i,j;
  realtype sum,sumi;

  if(sub_mean)set_mean();
  else set_copy();

  for(i=0;i<n;i++){
   sum=0.;
   sumi=0.;
   for(j=0;j+i<n;j++){
    sum+=arr[j]*arr[j+i];
    if(type==CORR_COMPLEX){  // (j)mal(j+i)*
     sum+=arri[j]*arri[j+i];
     sumi-=arri[j]*arr[j+i]-arr[j]*arri[j+i];
    }
   }
   tmp_arr[i]=sum/(n-i);
   if(type==CORR_COMPLEX)tmp_arri[i]=sumi/(n-i);
  }

  if(type==CORR_COMPLEX){

   for(i=0;i<n;i++){
    arr[i] =tmp_arr[i];
    arri[i]=tmp_arri[i];
   }
  }
  else{
   for(i=0;i<n;i++)arr[i]=tmp_arr[i];
  }

  delete [] tmp_arr;
  if(type==CORR_COMPLEX)delete [] tmp_arri;
  status=1;
}

void Correlation::insert_end(){
 // correcting possible error with the last value
 if(fini.stpi>0 && fini.stpi<n){
    //corr.insert_next(arr[corr.stpi-1]);
    fini.insert_next(ini[fini.stpi-1]);
 }
 if(type==CORR_COMPLEX){
   if(finii.stpi>0 && finii.stpi<n){
    //corri.insert_next(arri[corri.stpi-1]);
    finii.insert_next(inii[finii.stpi-1]);

   }
 }
 if(fini.stpi<n){
  printf("Correlation: warning: incomplete array set,\n"
         "padding %d values with zeroes.\n",n-fini.stpi);

  int i;
  for(i=fini.stpi;i<n;i++){
   //arr[i]=0.;
   ini[i]=0.;
   if(type==CORR_COMPLEX){
    //arri[i]=0.;
    inii[i]=0.;
   }
  }
 }
}


void Correlation::strait(){
 int i;
 if(sub_mean)set_mean();
 else set_copy();


 if(type==CORR_COMPLEX){ // a*(0) mal a(t)
  realtype x=arr[0];
  realtype y=arri[0];

  for(i=0;i<n;i++){
   arr[i] =arr[i]*x+arri[i]*y;
   arri[i]=arri[i]*x-arr[i]*y;
  }
 }
 else{
  realtype x=arr[0];
  for(i=0;i<n;i++)arr[i]=arr[i]*x;
 }
}


realtype Correlation::normalize(){
 realtype x;
 int i;
 if(type==CORR_COMPLEX){ // complex divide
  x=arr[0];
  //realtype y=arri[0];
  //realtype sq=x*x; //+y*y;
  for(i=0;i<n;i++){
   arr[i]/=x;
   arri[i]/=x;
  }
 }
 else{
  x=arr[0];
  for(i=0;i<n;i++)arr[i]=arr[i]/x;
 }
 return x;
}






double RunAverage::pnext(int dn){
  int i;
  double r;
  for(i=0;i<dn;i++){ // deleting the old data

    if(cur_n>=window){
      r=array[cur_n%window];
      s-=r;
      s2-=r*r;
      nsteps--;
    }
    array[cur_n%window]=cur; // adding new data
    cur_n++;
  }

  return Statistics2::pnext(dn);
}


int RunAverage::pushout(int dn){
  long k1, i;
  double r;

  if(dn>=window){
    if(cur_n<window)k1=cur_n;
    else k1=window;
    cur_n=0;
    clear();
    return k1;
  }

  k1=cur_n-window;
  if(k1<0)k1=0;
  long k2=cur_n-dn;
  if(k2<0)k2=0;

  cur_n-=dn;
  if(cur_n<0)cur_n=0;

  for(i=k1;i<k2;i++){
    r=array[i%window];
    s-=r;
    s2-=r*r;
    nsteps--;
    array[i%window]=array[(i+dn)%window];
  }
  return k2-k1;
}


int RunAverage::set_window(long w){
  if(w<=0)w=1;
  double *tmp=new double[w];
  if(!tmp){
    msg_error("RunAver: allocation error!\n"); 
    return 0;
  }

  long k1,k2,i;
  double r;

  if(w>window){
    k1=cur_n-window;
    k2=k1;
  }
  else{
    k1=cur_n-w;
    k2=cur_n-window;
  }
  if(k1<0)k1=0;
  if(k2<0)k2=0;

  if(w>window){
    cur_n=cur_n%window;
    k1=k1%window;
    k2=k2%window;
  }
 
  for(i=cur_n-1;i>=k1;i--){
    tmp[(i-k2)%w]=array[i%window];
  }
  for(;i>=k2;i--){
    r=array[i%window];
    s-=r;  // deleting old data
    s2-=r*r;
    nsteps--;
  }
  delete[] array;
    
  cur_n=(cur_n-k2)%w;
  window=w;
  array=tmp;

  return 1;
}
 

double RunAverage::rav() const {
  int i;
  double sum=0.;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  
  for(i=0;i<nn;i++){
    sum+=array[(cur_n-1-i)%window];
  }
  if(nn)return sum/nn;
  else return 0.;
}

double RunAverage::rav2() const {
  int i;
  double sum=0.;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  
  for(i=0;i<nn;i++){
    sum+=array[(cur_n-1-i)%window]*array[(cur_n-1-i)%window];
  }
  if(nn)return sum/nn;
  else return 0.;  
}


double RunAverage::rmin() const{
  int i;
  double sum=0.;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  if(nn==0)return 0.;
  double m=array[(cur_n-1)%window];
  for(i=1;i<nn;i++){
    double &v=array[(cur_n-1-i)%window];
    if(v<m)m=v;
  }
  return m;
}

double RunAverage::rmax() const {
  int i;
  double sum=0.;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  if(nn==0)return 0.;
  double m=array[(cur_n-1)%window];
  for(i=1;i<nn;i++){
    double &v=array[(cur_n-1-i)%window];
    if(v>m)m=v;
  }
  return m;
}

Statistics2& RunAverage::operator*=(double r){
  int i;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  
  for(i=0;i<nn;i++){
    array[(cur_n-1-i)%window]*=r;
  }
  
  return Statistics2::operator*=(r);
}


Statistics2& RunAverage::operator/=(double r){
  int i;
  int nn;
  if(cur_n<window)nn=cur_n;
  else nn=window;
  
  for(i=0;i<nn;i++){
    array[(cur_n-1-i)%window]/=r;
  }
  
  return Statistics2::operator/=(r);
}


  
  



