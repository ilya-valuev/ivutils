# ifndef _WINPLOT_H
# define _WINPLOT_H

# include <stdarg.h>
# include "Plot.h"

class WinPlot: public Plot{
 TForm *form;
 int hjust, vjust;
 int drback;
public:
 WinPlot(TForm *frm): Plot(){
  form=frm;
  hjust=0, vjust=0;
  drback;
 };

 void line(int x1, int y1, int x2, int y2){
  form->Canvas->MoveTo(x1,y1);
  form->Canvas->LineTo(x2,y2);
 }
 void setsolidbrush(int col){
  form->Canvas->Brush->Style=bsSolid;
  form->Canvas->Brush->Color=(TColor)col;
 }
 void bar(int x1, int y1, int x2, int y2){
  TRect r;
  r.Left=x1;
  r.Top=y1;
  r.Right=x2;
  r.Bottom=y2;
  form->Canvas->FillRect(r);
 }
 void setcolor(int col){
  form->Canvas->Brush->Color=(TColor)col;
  form->Canvas->Pen->Color=(TColor)col;
 }
 void rectangle(int x1, int y1, int x2, int y2){
  TRect r;
  r.Left=x1;
  r.Top=y1;
  r.Right=x2;
  r.Bottom=y2;
  form->Canvas->FrameRect(r);
 }

 void moveto(int x, int y, int draw=0){
  if(draw)form->Canvas->LineTo(x,y);
  else form->Canvas->MoveTo(x,y);
 }

 void moverel(int x, int y, int draw=0){
  if(draw)form->Canvas->LineTo(form->Canvas->PenPos.x+x,form->Canvas->PenPos.y+y);
  else form->Canvas->MoveTo(form->Canvas->PenPos.x+x,form->Canvas->PenPos.y+y);
 }

 int textwidth(char *str){
   return form->Canvas->TextWidth(AnsiString(str));
 }

 int textheight(char *str){
   return form->Canvas->TextHeight(AnsiString(str));
 }


 void TextJustifyLeftTop(){
   hjust=0;
   vjust=0;
 }

 void NoTextBackground(){
   drback=0;
 }

 void TextJustifyCenterTop(){
   hjust=1;
   vjust=0;
 }

 int grprintxy(int x,int y,char *format,...){
   if(drback)form->Canvas->Brush->Style=bsSolid;
   else form->Canvas->Brush->Style=bsClear;

   va_list ap;
   char string[250];
   int t;
   va_start(ap,format);
   t=vsprintf(string,format,ap);

   int h=0, v=0;
   if(hjust)h=-textwidth(string)/2;
   if(vjust)v=-textheight(string)/2;

   form->Canvas->TextOut(x+h,y+v,AnsiString(string));

   form->Canvas->Brush->Style=bsSolid;
   return t;
 }

 void TextJustyfyLeftCenter(){
   hjust=0;
   vjust=1;
 }

 void circle(int x, int y, int r){
   form->Canvas->Ellipse(x-r,y-r,x+r,y+r);
 }

};

# endif


