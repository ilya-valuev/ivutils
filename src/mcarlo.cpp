/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.8 $
 *   $Date: 2014/07/21 13:33:17 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/mcarlo.cpp,v 1.8 2014/07/21 13:33:17 valuev Exp $
 *
 *****************************************************************************/
#include<math.h>
#include<stdlib.h>
#include "mcarlo.h"
#include "common.h"


int MonteCarlo::test(realtype dE, realtype prefact){
 tested++;
 nstp++;
 int ac=0;
 if(dE<0 && prefact==1.)ac=1;
 else{
   realtype r1=prefact*exp(-dE/T);
   realtype r2=((double )rand())/RAND_MAX;
   if(r1<1){
     if(r1>=r2)ac=1;
   }
   else
     ac=1;
   //printf("+ %f ? %f\n",r1,r2);
 }

 if(ac){
  accepted++;
  accept++;
 }

 if(nstp%nav==0){
   aold=accept;
   accept=0;
   nstp=nav;
 }

 ratio=((double)(accept+aold))/((double)nstp);

 realtype kmax=1.1,kmax_inv=1./kmax;     // step auto-adjusting coefficient k
 if(nstp%nadj==0){
   need_adj=1;
 }
 if(ratio<0.5*kmax_inv)k=kmax_inv;
 else if(ratio>0.5*kmax)k=kmax;
 else k=ratio/0.5;
 
 return ac;
}

realtype MonteCarlo::adjust(){
 if(need_adj){
  need_adj=0;
  return k;
 }
 return 1.;
}

