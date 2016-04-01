# ifndef _PLOT_H
# define _PLOT_H


# include "common.h"

# define NO_DRAW 0xff8
/* if (back==NO_DRAW) don't draw background */


/* sides for arrows  */

# define LEFT 1
# define RIGHT 2
# define UP 1
# define DOWN 2

/* graphic line types */
# define SOLID 3
# define PUNCTIR 0
# define CIRCLES 1
# define CROSSES 2
# define GHISTOGRAM 4


/* line's clipping */
# define CUT_ON 1
# define CUT_OFF 0


/* types of variable's interdependence  */
# define FUNCTION 0
# define PARAM_FUNC 4
# define ARRAY 1
# define CHAIN 2
# define LINE 3
# define ONE_POINT 5




class Plot {
 int act,indx,indy,x,y,xs,xe,ys,ye;
 float x1,x2,y1,y2;
 float kx,ky;
 float xcur,ycur;
 int no_repeat;

 int lcolor, ltype,lcut;
 float lpar1,lpar2;


 int scrx(float x);
 int scry(float y);
 void check(float *x,float *y);
 void no_graph(){
   msg_error("Graphics not linked !\n");
 }

public:
 int pinit(int x_scr_start,int y_scr_start,int x_scr_end,int y_scr_end,
           float x1,float x2,float y1,float y2);
 Plot(int x_scr_start,int y_scr_start,int x_scr_end,int y_scr_end,
      float x1,float x2,float y1,float y2){
  pinit(x_scr_start,y_scr_start,x_scr_end, y_scr_end,x1,x2,y1,y2);
 

 }
 Plot(){
  act=0;
 }
 void pfield(int col,int back,int bord,int arrow_x,int arrow_y);
 void psetnames(char *xname,int xplace,char *yname,int yplace);
 void pgridx(float step,int hei);
 void pgridy(float step,int hei);
 void pscalex(float step,int hei,char *format);
 void pscaley(float step,int hei,char *format);
 void setplotstyle(int col,int type,...);
 void psetcolor(int col);
 int  pgetcolor(void);
 void pdrawplot(int type,...);

 void psetnorepeat(int s){
  no_repeat=s;
 }

 /* function for printing in graphics mode */
 # pragma argsused
 virtual int grprintxy(int x,int y,char *format,...){
  no_graph();
  return 0;
 }
 # pragma argsused
 virtual void line(int x1, int y1, int x2, int y2){
  no_graph();
 }
 # pragma argsused
 virtual void setsolidbrush(int col){
  no_graph();
 }

 # pragma argsused
 virtual void bar(int x1, int y1, int x2, int y2){
  no_graph();
 }

 # pragma argsused
 virtual void setcolor(int col){
  no_graph();
 }

 # pragma argsused
 virtual void rectangle(int x1, int y1, int x2, int y2){
  no_graph();
 }

 # pragma argsused
 virtual void moveto(int x, int y, int draw=0){
  no_graph();
 }

 # pragma argsused
 virtual void moverel(int x, int y, int draw=0){
  no_graph();
 }

 virtual void lineto(int x, int y){
  moveto(x,y,1);
 }

 virtual void linerel(int x, int y){
  moverel(x,y,1);
 }

 # pragma argsused
 virtual int textwidth(char *str){
  no_graph();
  return 0;
 }

 # pragma argsused
 virtual int textheight(char *str){
  no_graph();
  return 0;
 }

 virtual void TextJustifyLeftTop(){
   no_graph();
 }

 virtual void NoTextBackground(){
   no_graph();
 }

 virtual void TextJustifyCenterTop(){
   no_graph();
 }

 virtual void TextJustyfyLeftCenter(){
   no_graph();
 }

 # pragma argsused
 virtual void circle(int x, int y, int r){
   no_graph();
 }


};


# ifdef BGI

# include <graphics.h>
# include "grprint.h"

class BGIPlot: public Plot{
public:
 BGIPlot():Plot(){}
 BGIPlot(int x_scr_start,int y_scr_start,int x_scr_end,int y_scr_end,
         float x1,float x2,float y1,float y2):
    Plot(x_scr_start,y_scr_start,x_scr_end, y_scr_end,x1,x2,y1,y2){};

 void line(int x1, int y1, int x2, int y2){
  ::line(x1,y1,x2,y2);
 }
 void setsolidbrush(int col){
  setfillstyle(SOLID_FILL,col);
 }
 void bar(int x1, int y1, int x2, int y2){
  ::bar(x1,y1,x2,y2);
 }
 void setcolor(int col){
  ::setcolor(col);
 }
 void rectangle(int x1, int y1, int x2, int y2){
  ::rectangle(x1,y1,x2,y2);
 }

 void moveto(int x, int y, int draw=0){
  if(draw)::lineto(x,y);
  else ::moveto(x,y);
 }

 void moverel(int x, int y, int draw=0){
  if(draw)::linerel(x,y);
  else ::moverel(x,y);
 }

 int textwidth(char *str){
   return ::textwidth(str);
 }

 int textheight(char *str){
   return ::textheight(str);
 }

 void TextJustifyLeftTop(){
   settextjustify(LEFT_TEXT,TOP_TEXT);
 }

 void NoTextBackground(){
   setbackgroundmode(OFF);
 }

 void TextJustifyCenterTop(){
   settextjustify(CENTER_TEXT,TOP_TEXT);
 }

 int grprintxy(int x,int y,char *format,...){
   char str[200];
   va_list ap;
   va_start(ap,format);
   vsprintf(str,format,ap);
   return ::grprintxy(x,y,str);
 }

 void TextJustyfyLeftCenter(){
   settextjustify(LEFT_TEXT,CENTER_TEXT);
 }

 void circle(int x, int y, int r){
   ::circle(x,y,r);
 }

};


# endif






# endif