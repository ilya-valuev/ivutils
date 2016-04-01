# ifndef _ROOTPLOT_H
# define _ROOTPLOT_H


# include "matrix.h"
# include "plot.h"



class RootPlot{
 rootP Roots;
 int alloc;
public:
 int nroots;
 Plot *plt;

 float intx;
 float x1, x0, y1;
 float dx, dy;
 int th;
 float xav;
 float acc;
 int Ne;
 double Energy, Gap;
 int mparam;
 int colf,coln,colgr,colbk;

 RootPlot(Plot *p){
  plt=p;
  Roots=NULL;
  nroots=0;
  x0=-60;
  x1=30;
  y1=2;
  xav=0.001;
  dx=10.;
  dy=1.;
  th=3;
  alloc=0;
  Ne=0;
  mparam=0;
  colf=0x7,coln=0x1,colgr=0xfff,colbk=0x0f0;
 }
 ~RootPlot(){
  if(alloc)delete Roots;
 }
 void SetPointer(int n, rootP r){
  if(alloc)delete Roots;
  nroots=n;
  Roots=r;
  alloc=0;
 }

 int Read(char *file){
  if(alloc)delete Roots;
  nroots=0;
  alloc=0;
  Roots=NULL;
  if(ReadRoots(nroots, Roots,&Ne, file, -2)>=0)alloc =1;
  return alloc;
 }

 void SetInterval(float xs, float xe){
  if(x1<x0){
   x0=xe;
   x1=xs;
  }
  else{
   x0=xs;
   x1=xe;
  }
 }
 void AdjustXGrid(float xs, float xe){
  SetInterval(xs,xe);
  intx=(x1-x0)/20;
  float d=floor(log10(x1-x0));
  dx=2*(int)((x1-x0)*pow10(-d));
  dx=dx*pow10(d-1);
 }

 void SetXav(float x){
  xav=x;
  if(xav<0)xav=th/610.;
  acc=(x1-x0)*xav;
 }

 virtual void SetParams(){}

 void DisplayRoots(float left,float right,float av,int th,int filled);

 void Display(int sx1, int sy1, int sx2, int sy2){
  AdjustXGrid(x0,x1);
  SetXav(xav);

  plt->pinit(sx1,sy1,sx2,sy2,   x0,x1,0,y1);
  plt->pfield(colgr,colbk,0,1,1);

  /*setplotstyle(LIGHTRED,GHISTOGRAM);*/
  plt->setplotstyle(colf,SOLID);

  int filled=DetectParams(&Energy,&Gap,mparam,Ne,Roots);
  DisplayRoots(x0,x1,acc,th,filled);

  plt->pfield(colgr,NO_DRAW,0,1,1);
  plt->pgridx(dx,-3);
  plt->pgridy(dy,3);
  plt->pscalex(dx,4,"%2.2f");
 }
 # pragma argsused
 virtual void SetSolidPen(int thk){}
 # pragma argsused
 virtual void SetDottedPen(int thk){}

 void ShiftLeft(){
   x0-=intx;
   x1-=intx;
 }
 void ShiftRight(){
   x0+=intx;
   x1+=intx;
 }
 void ZoomIn(){
   float r=(x0+x1)/2,d=(x1-x0)/2;
   x1=r+d*1.1;
   x0=r-d*1.1;
 }
 void ZoomOut(){
   float r=(x0+x1)/2,d=(x1-x0)/2;
   x1=r+d/1.1;
   x0=r-d/1.1;
 }
 void Nplus(){
   y1+=1;
 }
 void Nminus(){
   y1-=1;
   if(y1<1.)y1=1.;
 }
 void AdjustMinMax();



};




# ifdef BGI

class BGIRootPlot: public RootPlot{
public:
 BGIRootPlot(Plot *p): RootPlot(p) {}
 void SetSolidPen(int thk){
  setlinestyle(SOLID_LINE,1,thk);
 }
 void SetDottedPen(int thk){
  setlinestyle(DOTTED_LINE,1,th);
 }
};


# endif

# if 0

class WinRootPlot: public RootPlot{
 TForm *frm;
public:
 WinRootPlot(TForm *f,Plot *p): RootPlot(p) { frm= f;}
 void SetSolidPen(int thk){
   frm->Canvas->Pen->Style=psSolid;
   frm->Canvas->Pen->Width=thk;
 }
 void SetDottedPen(int thk){
   frm->Canvas->Pen->Style=psDot;
   frm->Canvas->Pen->Width=thk;
 }
};

# endif





# endif