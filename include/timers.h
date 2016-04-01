/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : GridMD, ivutils, FDTD-II
 *
 *   $Revision: 1.2 $
 *   $Date: 2015/11/19 08:46:44 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/timers.h,v 1.2 2015/11/19 08:46:44 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/timers.h,v $
$Revision: 1.2 $
$Author: valuev $
$Date: 2015/11/19 08:46:44 $
*/
/*s****************************************************************************
 * $Log: timers.h,v $
 * Revision 1.2  2015/11/19 08:46:44  valuev
 * common base for AWP
 *
 * Revision 1.1  2009/04/29 13:45:37  valuev
 * added timers
 *
 *
*******************************************************************************/
/*r @file timers.h @brief Simple profiling tool 
*/ 

# ifndef _TIMERS_H
# define _TIMERS_H

# include "refobj.h"


class emTimer {
  refvector<emTimer> parentv;
  double time, t0;
  int stlev;
  
public:
  // the estimated processing time... is needed to turn off unnecessary timers
//  double etime;
//  string name;

  emTimer(int started=0, emTimer *parent_=NULL):parentv(0),time(0),t0(0)/*,etime(-1)*/,stlev(0){
    if(parent_)bind(parent_);
    if(started)start();
  }

  static double gettime(bool local=true);

  emTimer *bind(emTimer *p, int take_over=1);
 
  int start(double t0_=-1);

  int stop(int force=0, double tset_=-1);

  double update();
//  double elapsed() const {
//    if(!stlev)return gettime()-t0;
//    else return time;
//  }
  double *get_ptr(){
    return &time;
  }
  ~emTimer(){
    stop(1);
  }
};

# endif