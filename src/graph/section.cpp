# include "section.h"


int SectionCurve::AddPointVal(float x, float y, float val){
  if(AddPointRequest(x,y)){
    FillRequested(val);
  }
  return n;
}


void SectionCurve::to_integral(int num){
  int i;
  float dx,dy;
  for(i=np;i<n && i<np+num;i++){
    if(i==0)continue;
    dx=px[i]-px[i-1];
    dy=py[i]-py[i-1];
    sum+=lastval*sqrt(dx*dx+dy*dy);
  } 
}

float SectionCurve::Integral(){
  sum=0;
  int i;
  float dx,dy;
  for(i=1;i<np ;i++){
    dx=px[i]-px[i-1];
    dy=py[i]-py[i-1];
    sum+=0.5*(ampl[i]+ampl[i-1])*sqrt(dx*dx+dy*dy);
  } 
  return sum;
}


int SectionCurve::AddPointRequest(float x, float y){ 
  to_integral(n-np);

  for(;np<n;np++){
    pval.add(&lastval,1);
  }

  np=n;
  BoxedCurve::AddPoint(x,y);
  return n-np;
}


int SectionCurve::FillRequested(float val, int ind=-1){
  lastval=val;

  if(ind>=0 && np<n){
   pval.set(np-1,&val);
   to_integral(1);
   np++;
  }
  else{
    to_integral(n-np);
    for(;np<n;np++)pval.add(&val,1);
  }
  return n-np;
}


int SectionCurve::out(char *file){
   FILE *f=fopen(file,"wt");
   if(!f)return 0;

   fprintf(f,"#1-l/length 2-val/length 3-l 4-x 5-y 6-val 7-Dens(x) 8-Dens(y)\n");
   int i;
   float l=0,dx,dy,dl,densx,densy;

   for(i=0;i<np;i++){
     float a, b;
     if(fabs(length)<1e-32){
       a=0.;
       b=0;
     }
     else{
       a=l/length;
       b=ampl[i]/length;
     }
    
     if(i<np-1){
      dx=px[i+1]-px[i];
      dy=py[i+1]-py[i];
     }
     else{
       dx=0;
       dy=0;
     }
     dl=sqrt(dx*dx+dy*dy);

     if(fabs(dx)>1e-10){
       densx=ampl[i]*dl/dx;
     }
     else{
       densx=ampl[i]*dl*1e10;
     }

     if(fabs(dy)>1e-10){
       densy=ampl[i]*dl/dy;
     }
     else{
       densy=ampl[i]*dl*1e10;
     } 
     
     fprintf(f,"%e %e %e %e %e %e %e %e\n",
	     a,b,l,px[i],py[i],ampl[i],densx,densy);
     l+=dl;
   }
   fclose(f);
   return 1;
}


float SectionCurve::value_lnorm(float t){
  return value_l(t*length);
}


float SectionCurve::value_l(float l){
  if(np<2 || l<0 || l>length)return 0;
  
  int i;
  float suml=0,dx,dy,dl;
  for(i=1;i<np;i++){
    dx=px[i]-px[i-1];
    dy=py[i]-py[i-1];
    dl=sqrt(dx*dx+dy*dy);
    if(l<suml+dl){
      if(fabs(dl)<1e-32)return 0.5*(value(i)+value(i-1));
      else return value(i-1,(l-suml)/dl);
    }
    suml+=dl;
  }
  return value(np-1);
}


















