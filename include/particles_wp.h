/****************************************************************************
 *
 *   Copyright (c), Igor Morozov 2011        All Rights Reserved.
 *
 *   Author	: Igor Morozov, JIHT RAS, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *****************************************************************************/

#ifndef __PARTICLES_WP_H
#define __PARTICLES_WP_H

#include <vector>

#include "vector_3.h"

using namespace std;


// Class to store all particle data
class ParticlesWP {
public:
  enum {NONE = 0, SPLIT_WP = 0x1};

  unsigned flags;
  int ntot, ni, nup, ndown, ne, ne_real[2];
  vector<Vector_3> coord, vel, force;
  vector<double> width, pwidth;
  vector<double> mass, charge;
  Vector_3 *ri, *vi, *fi;
  Vector_3 *rdown, *vdown, *fdown, *rup, *vup, *fup;
  double *wdown, *wup, *pwdown, *pwup;
  Vector_2 *csplup, *cspldown;
  //int *nsplup, *nspldown;
  vector<int> nspl;       // number of splits for each electron
  vector<Vector_2> cspl;  // splitted WP coeffs
  int c_polar_mode;

  ParticlesWP(int ni_ = 0, int ndown_ = 0, int nup_ = 0, unsigned flags_ = NONE) {
    init(ni_, ndown_, nup_, flags_);
  }

  void init(int ni_ = 0, int ndown_ = 0, int nup_ = 0, unsigned flags_ = NONE) {
    flags = flags_;
    ni = ni_, ndown = ndown_, nup = nup_;
    ne = ndown + nup;
    ntot = ni + ne;

    if(ntot) {
      coord.resize(ntot);
      vel.resize(ntot);
      force.resize(ntot);
      width.resize(ntot);
      pwidth.resize(ntot);
      charge.resize(ntot);
      mass.resize(ntot);

      if(flags & SPLIT_WP) {
        cspl.resize(ne);
        nspl.resize(ne);
      }

      init_pointers();
    }

    c_polar_mode = 0;
  }

  void init_pointers() {
    if(ntot) {
      ri = &coord[0]; vi = &vel[0]; fi = &force[0];
    }
    else {
      ri = vi = fi = NULL;
    }

    if(ne) {
      rdown = ri + ni; vdown = vi + ni; fdown = fi + ni;
      rup = rdown + ndown; vup = vdown + ndown; fup = fdown + ndown;
      wdown = &width[0] + ni; wup = wdown + ndown;
      pwdown = &pwidth[0] + ni; pwup = pwdown + ndown;
      cspldown = &cspl[0]; csplup = cspldown + ndown;
    }
    else {
      rdown = vdown = fdown = NULL;
      rup = vup = fup = NULL;
      wdown = wup = NULL;
      pwdown = pwup = NULL;
      cspldown = csplup = NULL;
    }
    //nspldown = &nspl[0]; //nsplup = nspldown + ne_real[0];
  }

  ParticlesWP& operator=(ParticlesWP& p) {
    ni= p.ni; ndown = p.ndown; nup = p.nup; ne = p.ne; ntot = p.ntot; ne_real[0] = p.ne_real[0]; ne_real[1] = p.ne_real[1];
    flags = p.flags;
    coord = p.coord;
    vel = p.vel;
    force = p.force;
    width = p.width;
    pwidth = p.pwidth;
    charge = p.charge;
    mass = p.mass;

    if(flags & SPLIT_WP) {
      cspl = p.cspl;
      nspl = p.nspl;
    }

    init_pointers();

    return *this;
  }

  // Copy constructor
  ParticlesWP(ParticlesWP& p) {
    *this = p;
  };

  int copy_state(ParticlesWP& p, bool copy_ions = true) {
    if(ne != p.ne || ni != p.ni) return 1;

    if(copy_ions)
      for(int i=0; i<ni; i++) {
        ri[i] = p.ri[i];
        vi[i] = p.vi[i];
      }

    for(int i=0; i<ne; i++) {
      rdown[i] = p.rdown[i];
      vdown[i] = p.vdown[i];
      wdown[i] = p.wdown[i];
      pwdown[i] = p.pwdown[i];
    }

    if(flags & SPLIT_WP) {
      for(int i=0; i<ne; i++) cspl[i] = p.cspl[i];
      for(int i=0; i<ne_real[0]+ne_real[1]; i++) nspl[i] = p.nspl[i];
    }

    return 0;
  }

  double& get_par(int iwp, int ipar) {
    int i = iwp + ni;
    if(ipar < 3) return coord[i][ipar];
    if(ipar < 6) return vel[i][ipar-3];
    if(ipar == 6) return width[i];
    if(ipar == 7) return pwidth[i];
    return cspl[iwp][ipar-8];
  }
};


// Class for unit conversion
class UnitConverter {
public:
  int type;
  double dist, vel, acc, force, mass, time, eng, pres, charge;

  enum {IDENTITY = 1, TCPENGINE2ATOMIC, ATOMIC2TCPENGINE};

  UnitConverter(int type_ = 0) : type(0) {
    init(type_);
  }

  void init(int type_ = 0) {
    if(type_ == IDENTITY) {
       dist = vel = acc = force = mass = time = eng = pres = charge = 1.;
    }
    else if(type_ == TCPENGINE2ATOMIC) {
       dist = 1.8897259910357704;  // Angstrom -> Bohr radius (aB)
       vel = 0.19100733887787627;  // velocity units -> aB / atomic time units
       force = 0.01944689812235838;  // A*eV -> aB*Ehartr
       mass = 1.0072773476482406;  // m_proton -> atomic mass units
       eng = 0.03674930882684536;  // eV -> Hartree
       charge = 1.;
    }
    else if(type_ == ATOMIC2TCPENGINE) {
       dist = 1. / 1.8897259910357704;
       vel = 1. / 0.19100733887787627;
       force = 1. / 0.01944689812235838;
       mass = 1. / 1.0072773476482406;
       eng = 1. / 0.03674930882684536;
       charge = 1.;
    }
    else if(type_ > 0) type_ = 0;

    type = type_;
  }

  void clear() { type = 0; }
};

#endif
