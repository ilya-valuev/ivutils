# include "common.h"

int NewtonMaxIter= 1000;

int  NewtonSolve(float &xr,float acc,float func(float),float deriv(float),float xstart){
 
 float der,x,f;
 int sd;

 x=xstart;
 xr=xstart;
 der=deriv(x);
 f=func(x);
 if(fabs(f)<acc){
   xr=x;
   return 0;  // found
 }
 sd=SIGN(der);
 int iter=0;
 do{
   if(fabs(der)<1e-32)return -1;  // zero derivative
   //printf("%f %f %f\n",x,f,der);
   x=x-f/der;

   f=func(x);
   if(fabs(f)<acc){
    xr=x;
    return 0;  // found
   }

   der=deriv(x);
   if(sd!=SIGN(der)){
     return -2; // derivative changed sign
   }
   iter++;

 }while(iter<NewtonMaxIter);
 return -3; // too many iterations
}
 


