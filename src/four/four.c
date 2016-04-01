# include <math.h>
# include "four.h"
# include "common.h"

/* data[2*i] -- real part
   data[2*i+1] - imaginary part
   nn = 2^p,
   data must have 2*nn elements
   if( isign == 1)  direct FFT
   if( isign ==-1)  reverse FFT
   */

/*




*/
void four(fourtype *data, long nn, int isign){
 long i,istep,j,m,mmax,n;
 fourtype tempi,tempr;
 double theta,wi,wpi,wpr,wr,wtemp;

 n=2*nn;
 j=1;

 for(i=1;i<=n;i+=2){
  if(j>i){
   tempr=data[j-1];
   tempi=data[j];
   data[j-1]=data[i-1];
   data[j]=data[i];
   data[i-1]=tempr;
   data[i]=tempi;
  }

  m=n/2;
  while(m>=2 && j>m){
   j=j-m;
   m=m/2;
  }
  j=j+m;
 }

 mmax=2;
 while(n>mmax){
  istep=2*mmax;
  theta=(double)(6.28318530717959l/(isign*mmax));
  wpr=sin(0.5*theta);
  wpr=-2.*wpr*wpr;
  wpi=sin(theta);
  wr=1.;
  wi=0.;
  for(m=1;m<=mmax;m+=2){
   for(i=m;i<=n;i+=istep){
    j=i+mmax;
    tempr=((fourtype)wr)*data[j-1]-((fourtype)wi)*data[j];
    tempi=((fourtype)wr)*data[j]  +((fourtype)wi)*data[j-1];
    data[j-1]=data[i-1]-tempr;
    data[j]=data[i]-tempi;
    data[i-1]=data[i-1]+tempr;
    data[i]=data[i]+tempi;
   }
   wtemp=wr;
   wr=wr*wpr-wi*wpi+wr;
   wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
 }
}

void fold_it(fourtype *data, long nn){
  int i;
  for(i=1;i<=nn/2;i++){
    data[2*i]=(data[2*i]+data[2*nn-2*i])/nn;

    data[2*i+1]=(data[2*i+1]-data[2*nn-2*i+1])/nn;
    if(i!=nn/2){
      data[2*nn-2*i]=0;//data[2*i];
      data[2*nn-2*i+1]=0; //data[2*i+1];
    }
  }
  data[0]/=nn;
  data[1]/=nn;
}



fourtype *four_direct(fourtype *data, long nn, int isign, fourtype Ws, fourtype We, fourtype dW){

  fourtype w,cs,sn;
  long i,ic;
  int n=2*((int)((We-Ws)/dW)+5);
  fourtype *res=(fourtype *)malloc(n*sizeof(fourtype));
  if(!res)serror("four_direct: memory allocation error!\n");

  ic=0;
  for(w=Ws;w<=We;w+=dW){
    res[2*ic]=res[2*ic+1]=0;
    for(i=0;i<nn;i++){
      cs=(fourtype)cos(i*w);
      sn=(fourtype)sin(i*w);
      res[2*ic]+=cs*data[2*i]-sn*data[2*i+1];
      res[2*ic+1]+=sn*data[2*i]+cs*data[2*i+1];
    }
    res[2*ic]/=nn;
    res[2*ic+1]/=nn;
    ic++;
  }

  return res;
}





