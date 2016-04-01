/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.7 $
 *   $Date: 2015/10/16 18:47:00 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/mcarlo.h,v 1.7 2015/10/16 18:47:00 valuev Exp $
 *
 *****************************************************************************/
#ifndef mcarloH

# include "common.h"


class MonteCarlo {
protected:
 realtype T;
 int accept;
 int aold;
 int need_adj;
public:
 long tested;
 long accepted;
 int nstp;
 int nav; //e averaging interval, must be >> nadj
 int nadj; //e adjustment frequency
 realtype ratio;
 realtype k;

 void clear(){
  tested=accepted=0;
  accept=aold=0;
  nstp=0;
  ratio=0.;
  need_adj=0;
  k=1.;
 }
 MonteCarlo(){
  T=1.;
  nav=100;
  nadj=10;
  clear();
//# ifndef UNIX
//  randomize();
//# endif
 }

 realtype getT() const {
   return T;
 }
 
 void setT(const realtype &T_){
   T=T_;
 }

 long get_accepted() const{
   return accepted;
 }
 int test(realtype dE, realtype prefact=1.);
 realtype adjust();
};



#define mcarloH
//---------------------------------------------------------------------------
#endif
