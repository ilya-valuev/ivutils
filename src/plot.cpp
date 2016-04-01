
# include <conio.h>
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <stdarg.h>


# include "plot.h"
# include "common.h"


int Plot::pinit(int xs,int ys,int xe,int ye,float x1,float x2,float y1,float y2)
{
 act=0;
 xe=xe-8;
 ys=ys+8;
 if(xs>xe) swap(xs,xe,int );
 if(ys>ye) swap(ys,ye,int);
 if(x1>x2) swap(x1,x2,float);
 if(y1>y2) swap(y1,y2,float);
 if(x1==x2) { msg_error("pinit: X interval is zero.");return 0;}
 if(ys==ye) { msg_error("pinit: Y screen limit is zero.");return 0;}
 act=1;
 this->xs=xs;this->ys=ys;
 this->xe=xe;this->ye=ye;
 this->x1=x1;this->y1=y1;
 this->y2=y2;this->x2=x2;
 this->kx=(xe-xs)/(x2-x1);
 this->ky=(y2-y1)/(ye-ys);
 this->indy=this->indx=0;

 no_repeat=0;
 lcolor=0xffff;
 ltype=SOLID;
 lcut=CUT_OFF;
 act=1;
 return 1;

}



void Plot::pfield(col,back,bord,px,py)
int col,back,bord;
int px,py;
{
 int x,y,l,indx=0,indy=0;
 if(this->act){

  if(back!=NO_DRAW){
   setsolidbrush(back);
   //setfillstyle(SOLID_FILL,back);
   bar(this->xs,this->ys,this->xe,this->ye);
  }
  if(bord){
   setcolor(bord);
   rectangle(this->xs,this->ys,this->xe,this->ye);
  }

  x=this->xs-(this->kx)*(this->x1);
  y=this->ye+(this->y1)/this->ky;
  if(x>this->xe){ indx=1;x=this->xe;}
  if(x<this->xs){ indx=-1;x=this->xs;}
  if(y>this->ye){ indy=1;y=this->ye;}
  if(y<this->ys){ indy=-1;y=this->ys;}

  this->x=x;
  this->y=y;
  this->indx=indx;
  this->indy=indy;


  //if(!col)return;

  setcolor(col);
  l=y-4*abs(indx);
  if(l>=this->ys)line(x,this->ys,x,l);
  l=y+4*abs(indx);
  if(l<=this->ye)line(x,l,x,this->ye);
  l=x-4*abs(indy);
  if(l>=this->xs)line(this->xs,y,l,y);
  l=x+4*abs(indy);
  if(l<=this->xe)line(l,y,this->xe,y);

  if(indx*indy<=0){
   moveto(x,this->ys);
   linerel(0,-8);
   if(py){
    linerel(-2,4);
    moverel(4,0);
    linerel(-2,-4);
   }
   moveto(this->xe,y);
   linerel(8,0);
   if(px){
    linerel(-4,-2);
    moverel(0,4);
    linerel(4,-2);
   }
  }

 }
}



void Plot::psetnames(xname,xplace,yname,yplace)
char *xname,*yname;
int xplace,yplace;
{
 int w,h,x,y;
 w=textwidth(xname);
 h=textheight(xname);
 x=this->xe+6;
 if(xplace==UP) y=this->y-5-h;
 else y=this->y+4;
 if(this->indx>0)x=this->xs-w-2;
 //settextjustify(LEFT_TEXT,TOP_TEXT);
 TextJustifyLeftTop();
 grprintxy(x,y,xname);

 w=textwidth(yname);
 y=this->ys-10;
 if(yplace==LEFT)x=this->x-4-w;
 else x=this->x+12;
 if(this->indy<0)y=this->ye+2;
 grprintxy(x,y,yname);
}

void Plot::pgridx(float step,int hei)
{
 float steps;
 float xf=0,x=0;
 //int h;
 steps=step*this->kx;
 if(this->indx<=0){
  if(this->indx)xf=((int)(this->x1/step))*step;
  x=scrx(xf);/*this->x;*/
  for(;x<this->xe+1;x+=steps,xf+=step){
   if(x>=this->xs){
    moveto(x,this->y);
    linerel(0,hei);
   }
  }
 }
 xf=0.;
 if(this->indx>=0){
  if(this->indx)xf=((int)(this->x2/step))*step;
  x=scrx(xf); /*this->x;*/
  for(;x>=this->xs-1;x-=steps,xf-=step){
   if(x<=this->xe){
    moveto(x,this->y);
    linerel(0,hei);
   }
  }
 }
}

void Plot::pgridy(float step,int hei)
{
 float steps;
 float yf=0,y=0;
 steps=step/this->ky;
 if(this->indy>=0){
  if(this->indy)yf=((int)(this->y1/step))*step;
  y=scry(yf); /*this->y;*/
  for(;y>=this->ys-1;y-=steps,yf+=step){
   if(y<=this->ye){
    moveto(this->x,y);
    linerel(hei,0);
   }
  }
 }
 yf=0.;
 if(this->indy<=0){
  if(this->indy)yf=((int)(this->y2/step))*step;
  y=scry(yf); /*this->y;*/
  for(;y<this->ye+1;y+=steps,yf-=step){
   if(y>=this->ys){
    moveto(this->x+.5,y);
    linerel(hei,0);
   }
  }
 }
}


void Plot::pscalex(float step,int hei,char *format)
{
 float steps;
 float xf=0,x=0;
 int h,y;
 //setbackgroundmode(OFF);
 NoTextBackground();
 TextJustifyCenterTop();
 //settextjustify(CENTER_TEXT,TOP_TEXT);
 h=textheight("ILYA");
 steps=step*this->kx;
 y=this->y+hei;
 if(hei<0) y-=hei+2*h;
 if(this->indx<=0){
  if(this->indx)xf=((int)(this->x1/step))*step;
  x=scrx(xf);/*this->x;*/
  for(;x<this->xe+1;x+=steps,xf+=step){
   if(x>=this->xs)grprintxy(x,y,format,xf);
  }
 }
 xf=0.;
 if(this->indx>=0){
  if(this->indx)xf=((int)(this->x2/step))*step;
  x=scrx(xf); /*this->x;*/
  for(;x>=this->xs-1;x-=steps,xf-=step){
   if(x<=this->xe)grprintxy(x+.5,y,format,xf);
  }
 }
}

void Plot::pscaley(float step,int hei,char *format)
{
 float steps;
 float yf=0,y=0;
 int x;
 //setbackgroundmode(OFF);
 NoTextBackground();
 //settextjustify(LEFT_TEXT,CENTER_TEXT);
 TextJustyfyLeftCenter();

 steps=step/this->ky;
 x=this->x+2;
 //if(hei<0)
 x-=4+hei;
 if(this->indy>=0){
  if(this->indy)yf=((int)(this->y1/step))*step;
  y=scry(yf); /*this->y;*/
  for(;y>=this->ys-1;y-=steps,yf+=step){
   if(y<=this->ye)grprintxy(x,y,format,yf);
  }
 }
 yf=0.;
 if(this->indy<=0){
  if(this->indy)yf=((int)(this->y2/step))*step;
  y=scry(yf); /*this->y;*/
  for(;y<this->ye+1;y+=steps,yf-=step){
   if(y>=this->ys)grprintxy(x,y+.5,format,yf);
  }
 }
}





int Plot::pgetcolor(void){
 return lcolor;

}

void Plot::psetcolor(int c){
 lcolor=c;
}




void Plot::setplotstyle(int col,int type,...)
{
 va_list ap;
 char string[20];
 lcolor=col;
 lcut=CUT_OFF;
 ltype=type;
 if(type<3){
  va_start(ap,type);
  if(type==CROSSES){
   vsprintf(string,"%lf",ap);
   lpar1=atof(string)*this->kx;
   (int *)ap+=4;
   vsprintf(string,"%lf",ap);
   lpar2=atof(string)/this->ky;
  }
  else lpar1=(float)va_arg(ap,int);
 }
}

/*
float b(float x){
  getch();
  return(x);
} */

char string[30];

void Plot::check(float *x,float *y)
{
 if(*x < this->x1) *x=this->x1;
 if(*x > this->x2) *x=this->x2;
 if(*y < this->y1) *y=this->y1;
 if(*y > this->y2) *y=this->y2;
}

typedef float (*Func_ptr)(float);

void Plot::pdrawplot(int type,...)
{
 int j=0,l=0;
 char  *ty;
 char string[70];
 va_list ap;
 float dx,dy,q,t,w,x,xl,y,yl,x1,x2,step;
 /*,(*a)(float),(*b)(float);*/
 Func_ptr a,b;
 int ns,ne,nstep,size,on,k;
 int xsl,ysl,xs;
 //int ys;
 char *arr;
 va_start(ap,type);
 //ap=...;
 setcolor(lcolor);
 if(ltype==GHISTOGRAM){
  //setfillstyle(SOLID_FILL,lcolor);
  setsolidbrush(lcolor);
 } 

 switch(type){
  case ONE_POINT:vsprintf(string,"%lf %lf",ap);
		 sscanf(string,"%f %f",&x,&y);
		 if(this->act!=2){
		  xl=x;yl=y;
		  this->act=2;
		 }
		 else{
		  xl=this->xcur; yl=this->ycur;
		 }
		 break;
  case FUNCTION: vsprintf(string,"%p %lf %lf %lf",ap);
		 sscanf(string,"%p",&ty);
		 l=0;
		 while(string[l]!=' ')l++;
		 ns=sscanf(string+l,"%f %f %f",&x1,&x2,&step);
		 if(x1>x2) swap(x1,x2,float);
		 if(step<0) step=-step;
		 a=(Func_ptr)ty;
		 x=xl=x1;
		 yl=a(xl);
		 break;
  case PARAM_FUNC:
		 vsprintf(string,"%p %p %lf %lf %lf",ap);
		 sscanf(string,"%p %p",&ty,&b);
		 ns=sscanf(string+20,"%f %f %f",&x1,&x2,&step);
		 if(x1>x2) swap(x1,x2,float);
		 if(step<0) step=-step;
		 a=(Func_ptr)ty;
		 x=xl=a(x1);
		 yl=b(x1);
		 break;
  case CHAIN :
  case ARRAY :   vsprintf(string,"%p %d %d %d",ap);
		 sscanf(string,"%p %d %d %d",&arr,&ns,&ne,&nstep);
		 if(ns>ne) swap(ns,ne,int);
		 if(nstep<0) nstep=-nstep;
		 size=4;
		 if(type==CHAIN){
		  arr+=ns*size;
		  x=xl=ns;
		  y=yl=*((float *)arr);
		  arr+=nstep*size;
		  k=ns+nstep;
		 }
		 else{
		  arr+=2*ns*size;
		  x=xl=*((float *)arr);
		  arr+=size;
		  y=yl=*((float *)arr);
		  arr+=size+2*(nstep-1)*size;
		  k=ns+nstep;
		 }
		 break;
  case LINE  :   vsprintf(string,"%lf %lf %lf %lf",ap);
		 sscanf(string,"%f %f %f %f",&xl,&yl,&x,&y);
		 break;
 }
 on=1;


 moveto(q=scrx(xl),t=scry(yl));
 do{
  switch(type){
   case ONE_POINT: on=0;
		   break;
   case FUNCTION: x+=step;
		  if(x>x2){
		   x-=step;
		   on=0;
		  }
		  else y=a(x);
		  break;
   case PARAM_FUNC:
		  x1+=step;
		  if(x1>x2){
		   on=0;
		  }
		  else{
		   x=a(x1);
		   y=b(x1);
		  }
		  break;
   case ARRAY :   if(k>ne) on=0;
		  else{
		   x=*((float *)arr);
		   arr+=size;
		   y=*((float *)arr);
		   arr+=size+2*(nstep-1)*size;
		   k+=nstep;
		  }
		  break;
   case CHAIN :   if(k>ne) on=0;
		  else{
		   x=k;
		   y=*((float *)arr);
		   arr+=nstep*size;
		   k+=nstep;
		  }
		  break;
   case LINE :    on=0;
		  break;
  }
  if(no_repeat){
   if(x==xl)continue;
  }

  if(lcut){
   check(&x,&y);
   if(ltype==GHISTOGRAM) check(&xl,&yl);
   else{ xl=x;yl=y;}
  }

  switch(ltype){
   case SOLID:   lineto(scrx(x),scry(y));
		 break;
   case GHISTOGRAM:
		 if(yl==0.)break;
		 if(xl<x1)xl=x1;
		 if(xl>x2)xl=x2;
		 xsl=scrx(xl);
		 ysl=scry(yl);
		 xs=scrx(x);
		 //ys=scry(y);
		 dx=(xs-xsl)/2;
		 if(dx<1)dx=1;
		 bar(xsl-dx,this->y-1,xsl+dx,ysl);
		 /*bar(xsl+(xs-xsl)/2,this->y,xs,ys);*/
		 break;
   case CIRCLES: circle(scrx(xl),scry(yl),lpar1);
		 break;
   case CROSSES: moverel(lpar1 ,0);
		 linerel(-2*lpar1,0);
		 moverel(lpar1,lpar2);
		 linerel(0,-2*lpar2);
		 moveto(scrx(x),scry(y));
		 break;
   case PUNCTIR: w=sqrt((x-xl)*(x-xl)*(this->kx*this->kx)+(y-yl)*(y-yl)/(this->ky*this->ky));
		 if(w){
		  dx=lpar1*(x-xl)*this->kx/w;
		  dy=lpar1*(y-yl)/this->ky/w;
		  l=scrx(x);
		  /*w=scrx(xl); */
		  t=scry(yl);
		  while(q*( dx <0 ? -1 : 1) < l*(dx<0 ? -1 : 1)){
		   q+=dx;
		   t-=dy;
		   if(j) lineto(q,t);
		   else moveto(q,t);
		   j=1-j;
		  }
		 }
		 break;


  }
  xl=x;
  yl=y;
 }while(on);
 this->xcur=x;
 this->ycur=y;
}



int Plot::scrx(float x){
 return (this->xs+this->kx*(x-this->x1));
}

int Plot::scry(float c){
 return(this->ye-(c-this->y1)/this->ky);
}


# if 0

/*
float fx(float t){
 return(5*cos(t));
}

float fy(float t){
 return(-5*sin(t)/t);
} */


/*main(){
 int i,m=DETECT;
 int setplotstyle(int,int,...);
 char *buf;
 initgraph(&m,&m,"d:\\valuev\\turbo\\bgi");
 pinit(0,40,20,510,300,0.,50.,-1.,1.);
 pfield(0,WHITE,BLUE,RED,1,1);
 settextstyle(SMALL_FONT,HORIZ_DIR,USER_CHAR_SIZE);
 setusercharsize(1,1,1,1);
 psetnames("xss",DOWN,"yss",LEFT);
 /*pgridx(5.,p[n].ye-p[n].y-1);
 pgridy(2.,p[n].xe-p[n].x-1);
 pgridy(2.,-(p[n].x-p[n].xs-1));
 settextjustify(CENTER_TEXT,TOP_TEXT);
 pscalex(1.,-18,"%1.1lf");
 pgridx(.2,-2);
 pscaley(1.,2,"%1.1lf");   */
 setplotstyle(YELLOW,SOLID);
 for(l=0.;l<5.;l+=.01){
  xx=.9;
  setplotstyle(YELLOW,SOLID);
  pdrawplot(FUNCTION,ff,0.,50.,.5);
  delay(2);
  xx=.9;
  setplotstyle(BLUE,SOLID);
  pdrawplot(FUNCTION,ff,0.,50.,.5);
 }
} */

# endif


