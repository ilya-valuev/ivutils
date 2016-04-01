/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2012        All Rights Reserved.
 *
 *   Author  : Nikita Kazeev, Igor Morozov, Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project  : ivutils
 *
 *   $Revision: 1.5 $
 *   $Date: 2015/04/07 14:01:07 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/vts.h,v 1.5 2015/04/07 14:01:07 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/vts.h,v $
$Revision: 1.5 $
$Author: valuev $
$Date: 2015/04/07 14:01:07 $
*/
/*e****************************************************************************
 * $Log: vts.h,v $
 * Revision 1.5  2015/04/07 14:01:07  valuev
 * experiments with mode switching
 *
 * Revision 1.4  2014/06/23 16:13:57  morozov
 * Added old changes (on August 2013) related to the ion motion. Added box_hamiltonian files. Fixed compilation on VS2005 of cvector_3.h and box_hamiltonian.*.
 *
 * Revision 1.3  2012/08/30 09:04:07  morozov
 * forces_storage.h is replaced by common_def.h
 *
 * Revision 1.2  2012/06/27 14:09:14  morozov
 * VTS is optimized and moved into a separate class VTS_Control (to be implemented in LAMMPS)
 *
 * Revision 1.1  2012/06/21 09:15:33  valuev
 * added vts.h -- variable time stepping control
 *
 *
*******************************************************************************/
# ifndef VTS_H
# define VTS_H

/**\file vts.h \brief 
   \en Variable time stepping control for molecular dynamics.
   \ru —хема пременного временного шага дл€ молекул€рной динамики.
*/ 

#include "common_def.h"

enum VTS_METHOD {VTS_NONE = 0, VTS_SIMPLE, VTS_EPSILONS, VTS_PREDICTIVE};


///\en Provides variable time stepping control together with the logging at
///    equidistant time intervals
class VTS_Control {
protected:
  // internal VTS parameters
  enum VTS_METHOD method;
  enum STEP_ALGORITHM algorithm;
  double drift_max;
  double dt_natural;

public:
  // public VTS parameters to be changed directly
  double en_drift_goal; ///<\en desired average energy drift
  double en_drift_low;  ///<\en additional parameters for VTS_EPSILONS method: if
                        ///     en_drift_low < en_drift < en_drift_goal the timestep is unchnaged
  double dt_raise;      ///<\en maximal timestep increasing factor for a single step
  double dt_descent;    ///<\en timestep decreasing fator for VTS_SIMPLE and VTS_EPSILONS methods
  double exponent;      ///<\en exponent for VTS_PREDICTIVE method: dE ~ pow(dt,exp)
  double dt_low_limit;  ///<\en lowest allowed timestep, ignored if zero
  double dt_upper_limit; ///<\en highest allowed timestep, ignored if zero
  double dt_critical; ///<\en lowest step to autoswitch to ROLLBACK_ON_DE, ignored if zero


  ///\en sets the criteria for energy conservation at each step
  ///    ROLLBACK_ON_DE_DT the value of (E1-E0)/dt is preserved
  ///    ROLLBACK_ON_DE    the value of (E1-E0) is preserved
  enum ROLLBACK_MODE {ROLLBACK_ON_DE_DT, ROLLBACK_ON_DE} rollback_mode;

  // Counters used for diagnostocs
  int forwards;          ///<\en number of accepted steps
  int rollbacks;         ///<\en number of rollbacks
  int low_limit_reached; ///<\en number of approaching the lowest timestep

  ///\en Constructor. See \a init for explanation.
  VTS_Control(VTS_METHOD vts_method = VTS_NONE, STEP_ALGORITHM step_algorithm = STEP_MD) {
    init(vts_method, step_algorithm);
  }

  ///\en Sets the VTS method, basic stepping algorithm and initializes all 
  ///    public variables
  /// @param vts_method Timestep control method
  /// @param step_algorithm Global algorithm: MD or minimization
  void init(VTS_METHOD vts_method = VTS_NONE, STEP_ALGORITHM step_algorithm = STEP_MD);

  ///\en Sets the maximal value of energy change at which the rollback is requested 
  void set_max_drift(double en_drift_max = 0.);

  ///\en Returns \a true if a step made from t to t+dt is acceptable with the given energy drift \a en_drift.
  ///    In case of ACCEPTED step puts the suggested (possibly increased) NEXT time step
  ///    (to be used for stepping further from (t+dt) to (t+dt)+dt_new) into dt_new.
  ///    If time_to_update > 0, then dt_new is set to be lower or equal time_to_update. Used to match the 
  ///    time moments for log file updates.
  ///    In case of REJECTED step puts the suggested CORRECTED (reduced) time step (to be used for repeating a step
  ///    from t to t+dt_new) into dt_corrected. This value should be used for the next trial (revert to
  ///    time t, make step to t+dt_new, evaluate en_drift again).
  bool test_step(double dt, double en_drift, double* dt_new = 0, double time_to_update = 0.,
                 bool invalid_state = false);
};

# endif
