/****************************************************************************
 *
 *   Copyright (c), Igor Morozov 2011        All Rights Reserved.
 *
 *   Author	: Igor Morozov, JIHT RAS, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *****************************************************************************/

#include <wx/string.h>
#include <wx/tokenzr.h>

#include "state_io.h"
#include "logexc.h"


void HOOMDStateReader::GetUnits(const char* name) {
  wxString wxname(name);
  uc.clear();
  if( wxname == "atomic" )
    uc.init(UnitConverter::ATOMIC2TCPENGINE);
  else if( wxname == "tcpengine" )
    uc.init(UnitConverter::IDENTITY);
}


// Lookup the next node with the given name
wxXmlNode* HOOMDStateReader::FindNode(wxXmlNode* node, char* name, bool mandatory){
  while( node && node->GetName() != name )
    node = node->GetNext();

  if(mandatory && !node)
    LOGERR(-1,fmt("mandatory node %s isn't found", name),LINFO);

  return node;
}


// Lookup the next node with the given name
void HOOMDStateReader::read_particle_array(
  wxXmlNode* node, Vector_3* coord, double* width, Vector_2* cspl, int* nspl)
{
  GetUnits( node->GetPropVal("units", "") );
  if(!uc.type) LOGERR(-1,fmt("Unknown units in <%s>", node->GetName().c_str()),LINFO);

  int i, iwp = 0, counter[3] = {0, 0, 0};
  long idxe = -1, idxe_old = -1;

  if(cspl){
    for(i=0; i<ne; i++) nspl[i] = 0;  // clear split WP counters
    ne_real[0] = 0;
    ne_real[1] = 0;
  }

  wxStringTokenizer tkz( node->GetNodeContent() );
  for(i=0; tkz.HasMoreTokens(); i++) {
    if(i >= ntot)
      LOGERR(-1,fmt("too many particles in <%s> element", node->GetName().c_str()),LINFO);

    int itype = types[i];
    int idx = type_shift[itype] + counter[itype];  // index in the global particle array
    bool res = true;
    res = res   // read coordinates
      && tkz.GetNextToken().ToDouble(&coord[idx][0])
      && tkz.GetNextToken().ToDouble(&coord[idx][1])
      && tkz.GetNextToken().ToDouble(&coord[idx][2]);
    if(width) res = res && tkz.GetNextToken().ToDouble(&width[idx]);  // read width
    if(cspl) {  // if requested read WP spitting data
      Vector_2 c;
      res = res
        && tkz.GetNextToken().ToLong(&idxe)     // number of electron which WP belongs to
        && tkz.GetNextToken().ToDouble(&c[0])   // WP wavefunction coefficient (real)
        && tkz.GetNextToken().ToDouble(&c[1]);  // WP wavefunction coefficient (imag)

      if(idxe >=0 ) {  // ions are skipped due to idxe == -1
        if(idxe < idxe_old || idxe > idxe_old+1) LOGERR(-1,"Split WP is out of order",LINFO);
        
        nspl[idxe]++;
        cspl[iwp++] = c;
        if(idxe_old !=idxe)
          ne_real[itype-1]++;
        idxe_old = idxe;
      }
    }
    if(!res) LOGERR(-1,fmt("invalid numeric element for particle %d in <%s> element",
                    i, node->GetName().c_str()),LINFO);
    counter[itype]++;
  }

  //if(cspl) ne_real = idxe_old + 1;
}


// Open and analyse HOOMD input XML file
int HOOMDStateReader::Open(){
  // Open file
  if ( !doc.Load(m_filename) ) 
    return LOGERR(-1,"cannot open file",LINFO);

  // Check if the root node is correct
  root = doc.GetRoot();
  if (root->GetName() != "hoomd_xml") 
    return LOGERR(-1,"<hoomd_xml> is required",LINFO);

  return 0;
}


// Read particle types, coordinates and velocities
int HOOMDStateReader::ReadState(ParticlesWP& part){
  if(!root)
    return LOGERR(-1,"HOOMD file is not open",LINFO);

  // Get system type
  wxString aa("part_type"), bb("md");
  wxString ptype = root->GetPropVal("part_type", "md");
  if(ptype == "md") part_type = PTYPE_MD;
  else if(ptype == "wmpd") part_type = PTYPE_WPMD;
  else if(ptype == "wpmd_split") part_type = PTYPE_WPMD_SPLIT;
  else LOGERR(-1,"Invalid part_type in <hoomd_xml>",LINFO);

  // Retrieve the first (usually the sole) configuration
  wxXmlNode *conf =
    FindNode(root->GetChildren(), "configuration") -> GetChildren();

  // Get particle number and allocate structure
  wxStringTokenizer tkz(FindNode(conf, "type") -> GetNodeContent());
  ntot = tkz.CountTokens();

  // Create the table of particle types
  types.resize(ntot);
  int i, ni = 0, nedown = 0, neup = 0;
  for(i=0; tkz.HasMoreTokens(); i++) {
    wxString type = tkz.GetNextToken();
    if(type == "i") types[i] = 0, ni++;
    else if(type == "e") types[i] = 1, nedown++;
    else if(type == "f") types[i] = 2, neup++;
    else LOGERR(-1,fmt("unknown particle type '%s'", type.c_str()),LINFO);
  }
  ne = nedown + neup;

  // Allocate particle structure
  unsigned flags = 0;
  if(part_type >= PTYPE_WPMD_SPLIT) flags |= ParticlesWP::SPLIT_WP;
  part.init(ni, nedown, neup, flags);
  type_shift[0] = part.ri - &part.coord[0];
  type_shift[1] = part.rdown - &part.coord[0];
  type_shift[2] = part.rup - &part.coord[0];

  // Load positions, velocities and other properties depending on the particle type
  double *width = NULL, *pwidth = NULL;
  Vector_2* cspl = NULL;
  int* nspl = NULL;
  if(part_type >= PTYPE_WPMD) width = &part.width[0], pwidth = &part.pwidth[0];
  if(part_type >= PTYPE_WPMD_SPLIT) cspl = &part.cspl[0], nspl = &part.nspl[0];

  read_particle_array(FindNode(conf, "position"), &part.coord[0], width, cspl, nspl);
  read_particle_array(FindNode(conf, "velocity"), &part.vel[0], pwidth);

  if(part_type >= PTYPE_WPMD_SPLIT) {
    part.nspl.resize(ne_real[0]+ne_real[1]);
    part.ne_real[0] = ne_real[0];
    part.ne_real[1] = ne_real[1];
    // Read polar mode for complex coefficients
    long c_polar_mode;
    FindNode(conf, "position")->GetPropVal("c_polar_mode","0").ToLong(&c_polar_mode);
    part.c_polar_mode = c_polar_mode;
  }

  // Convert units
  for(i=0; i<ntot; i++) {
    part.coord[i]  *= uc.dist;
    part.width[i]  *= uc.dist;
    part.vel[i]    *= uc.vel;
    part.pwidth[i] *= uc.vel;
  }

  // Read masses
  wxXmlNode *node = FindNode(conf, "mass", false);
  if(node) {
    tkz.SetString( node->GetNodeContent() );
    int icount = 0;
    double value;
    for(i=0; tkz.HasMoreTokens(); i++) {
      wxString str = tkz.GetNextToken();
      if( !str.ToDouble(&value) )
        LOGERR(-1,fmt("invalid numeric element for particle %d in <mass> element",i),LINFO);
			if(icount >= ntot)
        LOGERR(-1,fmt("extra numeric element for particle %d in <mass> element",i),LINFO);
			else
        part.mass[icount++] = value*uc.mass;
    }
  } else {
    for(int i=0; i<ni; i++)
      part.mass[i] = 1836.1527556560675;  // assign hydrogen masses
    for(int i=0; i<ne; i++)
      part.mass[ni+i] = 1.;  // assign electron masses
  }

  // Read charges
  node = FindNode(conf, "charge", false);
  if(node) {
    tkz.SetString( node->GetNodeContent() );
    int icount = 0;
    double value;
    for(i=0; tkz.HasMoreTokens(); i++) {
      wxString str = tkz.GetNextToken();
      if( !str.ToDouble(&value) )
        LOGERR(-1,fmt("invalid numeric element for particle %d in <mass> element",i),LINFO);
      part.charge[icount++] = value*uc.charge;
    }
  } else {
    for(int i=0; i<ni; i++)
      part.charge[i] = 1.;  // assign ion charges
    for(int i=0; i<ne; i++)
      part.charge[ni+i] = -1.;  // assign electron charges
  }

  types.clear();

  return 0;
}


// Read box length
int HOOMDStateReader::ReadParameters(){
  if(!root)
    return LOGERR(-1,"HOOMD file is not open",LINFO);

  // Retrieve the first (usually the sole) configuration
  wxXmlNode* conf = FindNode(root->GetChildren(), "configuration");
  wxXmlNode* conf_child = conf->GetChildren();

  wxXmlNode* node = FindNode(conf_child, "box");
  GetUnits( node->GetPropVal("units", "") );
  if(!uc.type)  LOGERR(-1,fmt("Unknown units in <%s>", node->GetName().c_str()),LINFO);

  // Read box lengths
  if( !( 
    node->GetPropVal("Lx", "").ToDouble(&box[0]) &&
    node->GetPropVal("Ly", "").ToDouble(&box[1]) &&
    node->GetPropVal("Lz", "").ToDouble(&box[2]) ) )
    LOGERR(-1,"invalid numbers in <box>",LINFO);

  // Convert units
  for(int i=0; i<3; i++) box[i] *= uc.dist;

  // Read configuration time
  if( !conf->GetPropVal("time", "0.").ToDouble(&conf_time) )
    LOGERR(-1,"error reading time parameter in <configuration>",LINFO);

  // Read initial external energy
  if( !conf->GetPropVal("ext_energy", "0.").ToDouble(&conf_ext_energy) )
    LOGERR(-1,"error reading ext_energy parameter in <configuration>",LINFO);

  // Read initial step
  long step = 0;
  if( !conf->GetPropVal("time_step", "0").ToLong(&step) )
    LOGERR(-1,"error reading ext_energy parameter in <configuration>",LINFO);
  conf_time_step = (int) step;

  return 0;
}
