# include <math.h>
//# include <conio.h>
# include "vector_3.h"
# include "quater.h"



# if 0   //ndef _Windows

# include <graphics.h>

void ShowRotation(Basis_3& xyz_basis,void drawfunc(void),int swappages){
 Basis_3  rot(1,0,0,
              0,1,0,
              0,0,1);
 float alf;
 int page=0,ind,cont=0;
 char h;
 if(swappages){
  setactivepage(page);
  setvisualpage(page);
 }
 cleardevice();
 drawfunc();

 do{
  while(kbhit())h=getch();
  switch(h){
   case 'x':
    ind=1;
    break;
   case 'X':
    ind=-1;
    break;
   case 'y':
    ind=2;
    break;
   case 'Y':
    ind=-2;
    break;
   case 'z':
    ind=3;
    break;
   case 'Z':
    ind=-3;
    break;
   case ' ':
    cont=1-cont;
    h=0;
    break;
   default:
    ind=0;
  }
  if(ind){
   if(ind<0){
    ind=-ind;
    ind--;
    alf=-0.1; //-rot[ind].angle;
   }
   else { alf=0.1;--ind;} //rot[--ind].angle;
   xyz_basis=RotateBasis(xyz_basis,rot.vec[ind],alf);
   if(swappages)setactivepage(page=1-page);
   cleardevice();
   drawfunc();
   if(swappages)setvisualpage(page);
   ind=0;
   if(!cont)h=0;
  }
 }while(h!=27);
}

# endif