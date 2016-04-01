/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : GridMD, ivutils, FDTD-II
 *
 *   $Revision: 1.6 $
 *   $Date: 2014/07/18 14:04:37 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/timers.cpp,v 1.6 2014/07/18 14:04:37 morozov Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/timers.cpp,v $
$Revision: 1.6 $
$Author: morozov $
$Date: 2014/07/18 14:04:37 $
*/
/*s****************************************************************************
 * $Log: timers.cpp,v $
 * Revision 1.6  2014/07/18 14:04:37  morozov
 * Made compilable by icc 11.1
 *
 * Revision 1.5  2012/05/30 21:35:29  morozov
 * Fixed Unix version of emTimer::gettime
 *
 * Revision 1.4  2011/05/24 23:16:59  morozov
 * Using the real time clock instead of the CPU time
 *
 * Revision 1.3  2009/08/22 21:28:34  morozov
 * Corrected MPI version
 *
 * Revision 1.2  2009/08/22 12:58:00  morozov
 * Added '# ifdef USE_MPI' ...
 *
 * Revision 1.1  2009/04/29 13:45:37  valuev
 * added timers
 *
 *
*******************************************************************************/

# include <time.h>
# include "timers.h"

# ifdef USE_MPI
# include <mpi.h>
# endif


double emTimer::gettime(bool local){
#ifdef USE_MPI
  if(local)
    return MPI_Wtime();
  else{
    double sec;
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    if(irank==0)
      sec=MPI_Wtime();
    MPI_Bcast(&sec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    return sec;
  }
#else
#ifdef UNIX
  static timespec tspec;
  clock_gettime(CLOCK_REALTIME, &tspec);
  return tspec.tv_sec + 1.e-9 * tspec.tv_nsec;
#else
  return (double)clock()/(double)CLOCKS_PER_SEC;
#endif
#endif
}

emTimer *emTimer::bind(emTimer *p, int take_over){
  if(!take_over)return p;
  else{
    parentv.push_back(p);
    return this;
  }
}

int emTimer::start(double t0_){
  if(!stlev){
    t0 = t0_<0 ? gettime() : t0_;
    for(size_t i=0;i<parentv.size();i++)
      parentv[i]->start(t0);
  }
  return stlev++;
}

int emTimer::stop(int force, double tset_){
  if(!stlev)return 0; // not started
  if(force)stlev=0;
  else stlev--;
  if(!stlev){
    double tset = tset_<0 ? gettime() : tset_;
    time+=tset-t0;
    for(size_t i=0;i<parentv.size();i++)
      parentv[i]->stop(force, tset);
  }
  return stlev;
}

double emTimer::update(){
  if(stlev){
    double t1=gettime();
    time+=t1-t0;
    t0=t1;
  }
  return time;
}
