/****************************************************************************
 *
 *   Copyright (c), Igor Morozov 2011        All Rights Reserved.
 *
 *   Author	: Igor Morozov, JIHT RAS, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *****************************************************************************/

#ifndef __STATE_IO_H
#define __STATE_IO_H

#include <wx/xml/xml.h>
#include <wx/string.h>

#include "particles_wp.h"

class StateReader
{
protected:
  wxString m_filename;

public:
  StateReader(const char* fname = NULL) : m_filename(fname) {};
  virtual int Open() = 0;
  virtual int ReadState(ParticlesWP& part) = 0;
};


class HOOMDStateReader : public StateReader
{
protected:
  enum {PTYPE_MD = 0, PTYPE_WPMD, PTYPE_WPMD_SPLIT} part_type;
  int ntot, ne, ne_real[2];
  int type_shift[3];
  vector<int> types;
  wxXmlDocument doc;
  wxXmlNode *root;
  UnitConverter uc;

  void read_particle_array(wxXmlNode* node, Vector_3* coord,
    double* width = NULL, Vector_2* cspl = NULL, int* ispl = NULL);
  wxXmlNode* FindNode(wxXmlNode* node, char* name, bool mandatory = true);
  void GetUnits(const char* name);

public:
  Vector_3 box;
  double conf_time;
  double conf_ext_energy;
  int conf_time_step;

  HOOMDStateReader(const char* fname = NULL) : root(NULL), StateReader(fname), conf_ext_energy(0.), conf_time_step(0) {
    Open();
  };
  virtual int Open();
  virtual int ReadState(ParticlesWP& part);
  virtual int ReadParameters();
};

#endif
