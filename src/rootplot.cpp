
# include "rootplot.h"


void RootPlot::DisplayRoots(float left,float right,float av,int th,int filled){

 int i=0,n,k,passed=0;
 float curx;
 while(Roots[i].x<left && i<nroots)i++;
 if(i>filled){
  passed=1;
  plt->psetcolor(coln);
 }
 k=0;
 SetSolidPen(th);

 //mparam=1;
 int pparam;
 do{
  k++;
  if(i>=nroots)break;
  curx=left+k*av;
  n=0;

  if(i<=filled)pparam=1;
  else pparam=0;
  while(Roots[i].x<curx && i<nroots){n+=Roots[i++].n;};
  if(pparam){
   if(i<=filled+1){
    if(!mparam)pparam=0;
   }
  }

  if(n==0)continue;
  plt->pdrawplot(LINE,(float)curx-av/2,0.,(float)curx-av/2,(float)n);
  //if(!passed && i-1==filled)mparam=0;
  if(!passed && i>filled){
   passed=1;
   plt->psetcolor(coln);
   if(pparam){
    SetDottedPen(th);
    plt->pdrawplot(LINE,(float)curx-av/2,0.,(float)curx-av/2,(float)n);
    SetSolidPen(th);
   }
  }
 }while(curx<right);
 SetSolidPen(1);
}


void RootPlot::AdjustMinMax(){
  if(nroots==0)return;
  int i;
  float xmin=Roots[0].x, xmax=Roots[0].x;
  int nmax=Roots[0].n;

  for(i=1;i<nroots;i++){
    if(xmax<Roots[i].x)xmax=Roots[i].x;
    if(xmin>Roots[i].x)xmin=Roots[i].x;
    if(nmax<Roots[i].n)nmax=Roots[i].n;
  }
  if(xmin<0)x0=xmin*1.05;
  else x0=xmin/1.05;

  if(xmax>0)x1=xmax*1.05;
  else x1=xmax/1.05;

  y1=nmax+1;
 }



