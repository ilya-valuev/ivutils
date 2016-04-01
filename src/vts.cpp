/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2012        All Rights Reserved.
 *
 *   Author  : Nikita Kazeev, Igor Morozov, Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project  : ivutils
 *
 *   $Revision: 1.8 $
 *   $Date: 2015/04/07 14:01:07 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/vts.cpp,v 1.8 2015/04/07 14:01:07 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/vts.cpp,v $
$Revision: 1.8 $
$Author: valuev $
$Date: 2015/04/07 14:01:07 $
*/

/**\file vts.cpp \brief
   \en Variable time stepping control for molecular dynamics.
   \ru —хема пременного временного шага дл€ молекул€рной динамики.
*/

#include <algorithm>

#include <math.h>
#include <float.h>

#include "vts.h"

#ifdef _WIN32
#define isnan(x) _isnan(x)
#else
#define isnan(x) std::isnan(x)
#endif

using namespace std;


void VTS_Control::init(VTS_METHOD vts_method, STEP_ALGORITHM step_algorithm) {
  method = vts_method;
  algorithm = step_algorithm;
  rollbacks = forwards = low_limit_reached = 0;
  dt_natural = 0.;
  set_max_drift();

  // Default values for VTS parameters
  dt_low_limit = dt_upper_limit = 0.;
  dt_raise = (method == VTS_PREDICTIVE) ? 2.0 : 1.1;
  dt_descent = 0.5;
  exponent = 1.0;
  rollback_mode = ROLLBACK_ON_DE_DT;
	dt_critical = 0.;
}


void VTS_Control::set_max_drift(double en_drift_max) {
  drift_max = en_drift_max;
  // Default values for VTS parameters
  en_drift_goal = 0.2 * drift_max;
  en_drift_low = 0.1 * drift_max;
}


bool VTS_Control::test_step(double dt, double en_drift, double* dt_new, double time_to_update, bool invalid_state) {
  
  if(!method || !en_drift) {  // check main parameters
    if(dt_new) *dt_new = dt;
    return true;
  }

	if(dt_critical>0. && rollback_mode == ROLLBACK_ON_DE_DT && dt< dt_critical){  // switching mode
    rollback_mode = ROLLBACK_ON_DE;
    drift_max*= dt_critical;
		en_drift_goal*=dt_critical;
    en_drift_low*=dt_critical;
	}
	if(dt_critical>0. && rollback_mode == ROLLBACK_ON_DE && dt> dt_critical){  // switching mode
		rollback_mode = ROLLBACK_ON_DE_DT;
    drift_max/= dt_critical;
		en_drift_goal/=dt_critical;
    en_drift_low/=dt_critical;
	}


  double dE = fabs(en_drift);  // energy drift
  if(rollback_mode == ROLLBACK_ON_DE_DT) dE /= dt;  // taking dE/dt instead of dE
  int ferror = isnan(dE) || invalid_state;      // roll back on NaN values of energy

  bool accept = !( ferror
    || algorithm != STEP_MINIMIZE && dt > dt_low_limit && dE > drift_max
    || algorithm == STEP_MINIMIZE && en_drift > 0 );

  if(!dt_new) return accept;

  bool raise = accept;

  if(accept) {
    if (dt_natural) dt = dt_natural;
    // step is accepted
    forwards++;
    if (method == VTS_SIMPLE)
      dt *= dt_raise;
    else if (method == VTS_EPSILONS) {
      if (dE > en_drift_goal)
        dt *= dt_descent;
      else if (dE < en_drift_low)
        dt *= dt_raise;
    }
  } else {
    // step is rejected
    if (method == VTS_SIMPLE || method == VTS_EPSILONS || ferror) dt *= dt_descent;
    rollbacks++;
  }

  if (method == VTS_PREDICTIVE && !ferror) {
    double factor = min(pow(en_drift_goal/dE, 1./exponent), dt_raise);
    dt *= factor;
    raise = (factor > 1.);
  }

  if(raise) {  // timestep is increased
    // prevent from keeping too small timestep after adjusting to time_to_update
    if(dt < dt_natural) dt = dt_natural;
    else // Ensure hugh timestep limit
      if(dt_upper_limit && dt > dt_upper_limit) dt = dt_upper_limit;
  }
  else {  // timestep is decreased
    if(dt < dt_low_limit) {  // Ensure low timestep limit
    dt = dt_low_limit;
    low_limit_reached++;
    }
  }

  // Adjust timestep to the next update of the output files
  if(accept && time_to_update > 0) {
    dt_natural = dt;  // save original dt

    if (dt > time_to_update)
      dt = time_to_update;
    else if (dt > time_to_update/2)
      dt = time_to_update/2;
  }
  else
    dt_natural = 0.;  // indicates that the timestep was not adjusted

  *dt_new = dt;

  return accept;
}
