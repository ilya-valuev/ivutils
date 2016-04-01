# include <stdio.h>
# include <string.h>
# include <math.h>

# include "ncycle.h"



int CycleParam::init(char *spec,double start, double end, double step){
  strncpy(name,spec,250);
  p1=start;
  p2=end;
  dp=fabs(step);
  n=(int)ceil(fabs(p2-p1)/dp)+1;
  if(n<1)n=1;

  if(p2<p1)dp=-dp;

  i=0;
  nextp=NULL;

  // number of decimal digits required for oputput 

  int before=0,after=0;

  int l=(int)ceil(log10(fabs(p1)));
  if(l<=0)after=-l+1;
  else before=l;

  l=(int)ceil(log10(fabs(p2)));
  if(l<=0)after=fmax(-l+1,after);
  else before=fmax(l,before);

  l=(int)ceil(log10(fabs(dp)));
  if(l<=0)after=fmax(-l+1,after);
  else before=fmax(l,before);


  if(after>0)after+=2; // precision is 2 digits
  //printf("before: %d, after: %d\n",before,after);

  sprintf(digfmt,"%%0%d.%df",before,after);
  return n;
}




                   
int CycleParam::init(char *spec,char *range,char *delim=","){
  float a1,a2,s;
  int r=scan_range(range,&a1,&a2,&s,delim);
  if(r==3){
    return init(spec,a1,a2,s);
  }
  else if(r==1){
    return init(spec,a1,a1,1.);
  }
  else {
    msg_error("CycleParam: Invalid initialization: %s %s\n",spec,range);
  }
  return 0;
}





int Cycle::init(int n){
  maxnp=n;
  p= new CycleParamP[n];
  if(!p)return 0;
  int i;
  for(i=0;i<n;i++)p[i]=NULL;
  np=0;
  return 1;
}


int Cycle::Next(){
  if(np<1)return 0;

  int lm=1,sum=0;
  int i;
  (p[0]->i)++;
  for(i=0;i<np;i++){
   if(p[i]->i>=p[i]->n && i!=np-1){
    p[i]->i=0;
    (p[i+1]->i)++;
   }
   sum+=lm*p[i]->i;
   lm*=p[i]->n;
  }
  if(p[np-1]->i<p[np-1]->n)return sum+1;
  else return 0;
}


int Cycle::Current(){
  if(np<1)return 0;

  int lm=1,sum=0;
  int i;
  
  for(i=0;i<np;i++){
   sum+=lm*p[i]->i;
   lm*=p[i]->n;
  }
  return sum+1;
} 


int Cycle::Total(int lev){
  if(np<1)return 0;
  if(lev<0){
    int i,sum=0;
    for(i=0;i<np;i++)sum+=p[i]->n;
    return sum;
  }
  if(lev>=np)return 0;
  return p[lev]->n;
}


char *Cycle::SpecString(){
  char str[200];
  int i;
  
  for(i=np-1;i>=0;i--){
    strcpy(str,"");
    strcat(str,p[i]->getname());
    strcat(str,p[i]->format());
    if(i!=0)strcat(str,"_"); // delimiter
    sprintf(buffer,str,p[i]->value());
  }
  return buffer;
}






