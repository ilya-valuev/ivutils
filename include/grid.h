/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.6 $
 *   $Date: 2012/06/29 10:50:12 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/grid.h,v 1.6 2012/06/29 10:50:12 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/grid.h,v $
$Revision: 1.6 $
$Author: valuev $
$Date: 2012/06/29 10:50:12 $
*/
/*s****************************************************************************
 * $Log: grid.h,v $
 * Revision 1.6  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.45  2012/03/22 07:31:20  valuev
 * updated increments for iterators, prepared parallel/universal near-to-far (not working)
 *
 * Revision 1.44  2012/03/21 17:09:37  lesha
 * documentation
 *
 * Revision 1.43  2011/09/23 13:52:53  valuev
 * compiled with gcc 4.4
 *
 * Revision 1.42  2011/09/14 18:37:41  lesha
 * ch_drc is added. some modification with read table
 *
 * Revision 1.41  2011/05/05 22:53:18  lesha
 * working axis rules are changed
 *
 * Revision 1.40  2011/04/28 03:11:15  lesha
 * confine_memory modification etc. (transport)
 *
 * Revision 1.39  2011/04/27 16:56:42  lesha
 * fixing allocation memory problems in transport
 *
 * Revision 1.38  2011/03/18 20:29:15  lesha
 * make transport code compilable
 *
 * Revision 1.37  2011/03/18 17:51:00  biaks
 * Same errors fixed.
 *
 * Revision 1.36  2011/02/14 20:32:18  lesha
 * confinement is added to NonUniformGrid
 *
 * Revision 1.35  2011/02/12 23:49:28  lesha
 * dim is added to NonUniformGrid
 *
 * Revision 1.34  2011/02/01 18:59:27  lesha
 * nonuniform grid is included
 *
 * Revision 1.33  2010/10/24 23:31:48  valuev
 * added parallel VTK writer
 *
 * Revision 1.32  2010/06/06 11:15:30  lesha
 * const_iterator is added to vector_set
 *
 * Revision 1.31  2010/04/08 17:49:26  lesha
 * some comments are added
 *
 * Revision 1.30  2009/03/23 20:57:11  valuev
 * corrected interleaved loop
 *
 * Revision 1.29  2009/03/20 10:07:16  valuev
 * global memory interleave, updated dipoles
 *
 * Revision 1.28  2009/03/18 10:33:28  valuev
 * added support for interleaved grids
 *
 * Revision 1.27  2008/04/23 02:08:04  lesha
 * test_local is added
 *
 * Revision 1.26  2008/04/08 23:37:39  lesha
 * emFluxSpectra is added
 *
 * Revision 1.25  2008/01/07 01:10:08  lesha
 * make gcc compilable
 *
 * Revision 1.24  2008/01/06 03:25:31  lesha
 * Implementation part is moved from *.h to *.cpp
 *
 * Revision 1.23  2007/11/13 14:41:06  valuev
 * Added randomization for degenerated case in simplex;
 * Corrected some grammar
 *
 * Revision 1.22  2007/10/12 10:34:12  lesha
 * get_dx is added
 *
 * Revision 1.21  2007/08/13 21:42:49  lesha
 * Oblique correction is added to emDispersiveOffdiag (test version): all_nonlocal flag is added
 *
 * Revision 1.20  2007/06/27 23:25:13  lesha
 * new detectors parallelization
 *
 * Revision 1.19  2007/05/30 17:17:04  valuev
 * New emInterpolation implementation
 *
 * Revision 1.18  2007/05/29 22:32:05  lesha
 * valtype -> vec_type (for g++ compibility)
 *
 * Revision 1.17  2007/05/25 13:01:38  valuev
 * oblique incidence (test version)
 *
 * Revision 1.16  2007/04/19 22:17:14  valuev
 * modified normal detection algorithm
 *
 * Revision 1.15  2007/03/20 17:30:07  lesha
 * include common.h
 *
 * Revision 1.14  2007/03/13 06:20:35  lesha
 * function accdiv is added
 *
 * Revision 1.13  2007/03/11 00:48:12  lesha
 * VEC_ZERO2 = 1e-8
 *
 * Revision 1.12  2007/03/08 11:43:58  lesha
 * one error is corrected
 *
 * Revision 1.11  2007/03/07 23:29:10  lesha
 * Some errors are corrected
 *
 * Revision 1.10  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.9  2007/01/26 13:34:57  valuev
 * added transfer buffers
 *
 * Revision 1.8  2007/01/26 09:44:14  valuev
 * added external border, transfers are moved back to analyze_contours
 *
 * Revision 1.7  2007/01/24 09:14:45  valuev
 * corrected memory range
 *
 * Revision 1.6  2007/01/23 23:05:04  valuev
 * added SetMemoryRange
 *
 * Revision 1.5  2007/01/22 13:50:20  valuev
 * Added used/unused grid for Yee block, added channels for RecordAnalyzer
 *
 * Revision 1.4  2006/12/05 22:52:34  valuev
 * Fixed get_group_count, added zero length group support to group_pack
 *
 * Revision 1.3  2006/11/29 18:05:05  valuev
 * made the code compilable with g++
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
# ifndef _GRID_H
# define _GRID_H

#include <algorithm>
#include "region.h"

const int global_il=1;

/// uniform grid in space
template <class value_tt, size_t interleave=1>
class UniformGrid {
protected:
  value_tt *ptr;
  Vector_3 dx;
  Vector_3 pref;
  Vector_3 ds;
  Box b;
  /// interpolation settings: default
  int ngr[3];
  int nst[3];
  int nen[3];
  /// memory settings
  int mgr[3];
  int mst[3];
  int men[3];
  int SZ2; //< mgr[1]*mgr[2], for unpack_ind
  int GSZ;
  int nintz; // z-size of the interleaved data
  /// iterator settings
  int igr[3];
  int ist[3];
  int ien[3];

  int drc; // bit flag responsible for axis directions
  Vector_3 ch_drc(const Vector_3 &pos)const{
    Vector_3 p1=b.get_p1(),p2=b.get_p2(),p;
    for(int i=0;i<3;i++){
      if(drc&(1<<i))
        p[i]=p2[i]-(pos[i]-p1[i]);
      else
        p[i]=pos[i];
    }
    return p;
  }
  
  void set_range(int *st, int *en, int *gr, const int *start, const int *end);
public:
  typedef value_tt value_t;
  /// Iterator going through all points in the grid.
  /// Index change order is Z(2), Y(1), X(0), index goes from small to large
  class iterator{
  protected:
    const UniformGrid *parent;
    int ind[3];  
    int ie;
    /// grid constructor, constructs sequence start 
    iterator(const UniformGrid *sparent, int end_=0):parent(sparent),ie(end_){
      for(int i=0;i<3;i++)ind[i]=parent->ist[i];
    }
  public:
    typedef int difference_t;
    friend class UniformGrid;
    /// default constructor constructs sequence end
    iterator():parent(NULL),ie(1){
      for(int i=0;i<3;i++)ind[i]=0;
    }
    /// copy constructor 
    iterator(const iterator& other):parent(other.parent),ie(other.ie){
      for(int i=0;i<3;i++){
        ind[i]=other.ind[i];
      }
    }
    /// iterator difference
    int operator-(const iterator& other) const;

    Vector_3 operator*() const {
      return parent->Position(ind[0],ind[1],ind[2]);
    }

    Vector_3 ds() const{
      return parent->ds;
    }
    /// iterates contours in the same order as in explicit update
    /// index change order:  iz, iy, ix, dir (0,1,2), ftype (E=0,H=1)
    iterator& operator++();
     
    /// positive increment to an iterator
    iterator& operator+=(int incr);

    /// iterates contours in the same order as in explicit update
    /// index change order:  iz, iy, ix, dir (0,1,2), ftype (E=0,H=1)
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }

    /// complete test
    bool operator!=(const iterator &other) const{
      if(ie!=other.ie)
        return true;
      if(ie==1)
        return false; // ends always equal
      for(int i=0;i<3;i++){
        if(ind[i]!=other.ind[i])
          return true;
      }
      return false;
    }

    iterator get_grid_it(int gind){
      iterator o=*this;
      parent->unpack_ind(gind,o.ind[0],o.ind[1],o.ind[2]);
      return o;
    }

    iterator shift(int *sh){
      iterator o=*this;
      for(int i=0;i<3;i++){
        o.ind[i]+=sh[i];
      }
      return o;
    }

    int get_size() const {
      return get_size(0)*get_size(1)*get_size(2);
    }

    int get_size(int i) const {
      return parent->ien[i]-parent->ist[i]+1;
    }
    int get_end(int i) const {
      return parent->ien[i];
    }

    int get_start(int i) const {
      return parent->ist[i];
    }

    int get_ind(int i) const {
      return ind[i];
    }
  };
  typedef iterator const_iterator;

  iterator begin() const {
    return iterator(this);
  }

  iterator end() const {
    return iterator(this,1);
  }
  
  UniformGrid():ptr(NULL),drc(0){}

  // v1, v2 - противоположные вершины куба, ограничивающего сетку, 
  // sz - массив из трех чисел, указующих на количество шагов сетки по трем направлениям
  UniformGrid(const Vector_3 &v1, const Vector_3&v2, const int dir, const int *sz, const int *start=NULL, const int *end=NULL):ptr(NULL){
    init(v1,v2,dir,sz,start,end);
  }

  UniformGrid(const Vector_3 &v1, const Vector_3&v2, const int *sz, const int *start=NULL, const int *end=NULL):ptr(NULL){
    init(v1,v2,0,sz,start,end);
  }

  size_t get_interleave() const {
    return interleave;
  }

  void init(const Vector_3 &v1, const Vector_3&v2, const int dir, const int *sz, const int *start=NULL, const int *end=NULL);

  /// sets limits for iterator
  void SetIteratorRange(const int *start=NULL, const int *end=NULL){
    set_range(ist,ien,igr,start,end);
  }

  /// sets limits for iterator
  void SetInterpolationRange(const int *start=NULL, const int *end=NULL){
    set_range(nst,nen,ngr,start,end);
  }

  void GetInterpolationRange(int *gr){
    for(int i=0;i<3;i++)
      gr[i]=ngr[i];
  }

  /// set memory block ranges to be addressed
  /// if interleave is local and nonunity, memory size in direction 2 must be a multiple of the interleave 
  void SetMemoryRange(const int *start=NULL, const int *end=NULL){
    set_range(mst,men,mgr,start,end);
    SZ2=mgr[1]*mgr[2];
    GSZ=mgr[0]*mgr[1]*mgr[2];
    if(global_il){ // interleave the whole grid
      //int aux_cells=2;
      // maximal plane size
      //int szm=max(mgr[0]*mgr[1],SZ2);
      //szm=max(szm,mgr[0]*mgr[2]);
      //GSZ+=aux_cells*szm;
      int rest=GSZ%interleave;
      if(rest)
        GSZ+=interleave-rest;
      nintz=GSZ/interleave;
    }
    else{  // interleave the z-direction only 
      nintz=mgr[2]/interleave;
      if(mgr[2]%interleave)  // this will cause divide by zero error in interpolation to indicate the wrong setting
        nintz=0; 
    }
  }

  /// gets active array size in memory
  size_t msize(int dir=-1) const {
    return dir>=0 ? (size_t)mgr[dir] : (size_t)mgr[0]*mgr[1]*mgr[2];
  }

  /// gets full size in memory (including auxilliary cells)
  size_t size() const {
    return GSZ;
  }

  Vector_3 get_dx() const {
    return dx;
  }

  Vector_3 get_ds() const {
    return ds;
  }

  Vector_3 get_total_surface(){
    Vector_3 total_ds=ds;
    for(int i=0;i<3;i++)
      if(!ds[i])total_ds*=ngr[i];
    return total_ds;
  }

  /// @returns the position of a given grid point
  Vector_3 Position(int ix, int iy, int iz) const {
    return pref+Vector_3(ix*dx[0],iy*dx[1],iz*dx[2]);
  }

  inline int pack_base_ind(int ix, int iy, int iz) const {
    return ((ix-mst[0])*mgr[1]+iy-mst[1])*mgr[2]+iz-mst[2];
  }


  inline int pack_base_ind(int ix, int iy) const {
    return ((ix-mst[0])*mgr[1]+iy-mst[1])*mgr[2];
  }

  inline int pack_z_ind(int iz) const {
    if(interleave==1 || global_il)
      return iz-mst[2];
    else{      
      iz-=mst[2];
      return iz/nintz+interleave*(iz%nintz);
    }
  }

  inline int pack_z0_ind(int iz) const {
    if(interleave==1 || global_il)
      return iz;
    else{      
      return iz/nintz+interleave*(iz%nintz);
    }
  }


  inline int pack_ind(int ix, int iy, int iz) const {
    if(global_il && interleave!=1){
      int ind=pack_base_ind(ix,iy,iz);
      return ind/nintz+interleave*(ind%nintz);
    }
    else
      return pack_base_ind(ix,iy)+pack_z_ind(iz);
  }

  inline void unpack_ind(int ind, int &ix, int &iy, int &iz) const {
    if(global_il && interleave!=1)
      ind=(ind%interleave)*nintz+ind/interleave;
    ix=ind/SZ2+mst[0];
    ind%=SZ2;
    iy=ind/mgr[2]+mst[1];
    if(interleave==1 || global_il)
      iz=ind%mgr[2]+mst[2];
    else{
      iz=ind%mgr[2];
      iz=(iz%interleave)*nintz+iz/interleave+mst[2];
    }
  }

  inline int unpack_z_ind(int iz) const {
    if(interleave!=1 && !global_il)
      return (iz%interleave)*nintz+iz/interleave+mst[2];
    else
      return iz+mst[2];
  }

  inline int incr_z0_ind(int iz, int sh) const {
    if(interleave!=1){
      iz= (iz%interleave)*nintz+iz/interleave; // unpack
      iz+=sh; // shift
      if(global_il){
        //if(iz>=GSZ)
          //return -1;
        iz=iz%GSZ;
      }
      else
        iz=iz%mgr[2];
      return iz/nintz+interleave*(iz%nintz); // pack
    }
    else
      if(global_il)
        return (iz+sh)%GSZ;
      else 
        return (iz+sh)%mgr[2];
  }

  /// gets the loop ranges for two loops: [0, range1) and [range1, nintz)
  /// which must be performed for each interleave bank for
  /// convolution with given positive shift
  /// returns: bank section point: range1 
  ///          bank shift for the first loop (the other one is dbank0+1),
  ///          shifts for the first and the second loop
  void get_contiguous_ranges(int shift, int &range1, int &dbank0, int &shift0, int &shift1) const {
    dbank0=shift/nintz;
    shift0=shift%nintz;
    range1=nintz-shift0;
    shift1=-range1;
  }

  int get_range1(int shift) const {
    return nintz-shift%nintz;
  }

  /// validity check for index
  bool check_ind(int ix, int iy, int iz) const {
    if(ix<mst[0] || ix>men[0])return false;
    if(iy<mst[1] || iy>men[1])return false;
    if(iz<mst[2] || iz>men[2])return false;
    return true;
  }

  /// check for interpolation index
  bool check_interpolation_ind(int ix, int iy, int iz) const {
    if(ix<nst[0] || ix>nen[0])return false;
    if(iy<nst[1] || iy>nen[1])return false;
    if(iz<nst[2] || iz>nen[2])return false;
    return true;
  }

  value_tt &operator()(int ix, int iy, int iz) const {
    return ptr[pack_ind(ix,iy,iz)];
  }

  value_tt operator()(const Vector_3 &place) const{
    return Interpolate(place);
  }

  /// distributes the value being added  between nearest grid points 
  /// @returns the number of grid points affected
  int AddValue(const Vector_3 &place, value_tt value, int force_external=0);

  int SetPtr(value_tt *sptr){
    ptr=sptr;
    return 1;
  }

  /// Gets interpolation coefficients
  /// the arrays must be at least 8 elements long
  /// if nonlocal is not NULL,  
  /// the negative indicies in int_ind correspond to absent grid points
  /// the space location of kth absent point is pushed_back to nonlocal[-int_ind[k]-1]
  int GetCoeff(const Vector_3 &place, int *int_ind, vec_type *values, int force_external=0, vector<Vector_3> *nonlocal=NULL, bool all_nonlocal=false) const;

#if 0
  /// gets interpolation indicies and coefficients
  /// indicies are in the relative form (shifts) from major_index
  /// @return major_index or -1 if out of region
  int GetShifts(const Vector_3& place, vector<int> &shifts, vector<vec_type> &coeffs, int force_external=0) const{
    int ind[8];
    vec_type val[8];
    int nv=GetCoeff(place,ind,val,force_external);
    if(nv<=0)return -1;
    shifts.push_back(nv);
    int i, m_ind=ind[0];
    for(i=0;i<nv;i++){
      ind[i]-=m_ind;
      shifts.push_back(ind[i]);
      coeffs.push_back(val[i]);
    }
    return m_ind;
  }
 #endif

  int test_local(const Vector_3 &place) const;

  value_tt Interpolate(const Vector_3 &place, int force_external=0) const;
};

template <class value_tt, size_t interleave=1>
class NonUniformGrid: public UniformGrid<value_tt,interleave>{

  int dim_type[3];

  using UniformGrid<value_tt,interleave>::ptr;
  using UniformGrid<value_tt,interleave>::dx;
  using UniformGrid<value_tt,interleave>::pref;
  using UniformGrid<value_tt,interleave>::ds;
  using UniformGrid<value_tt,interleave>::b;
  /// interpolation settings: default
  using UniformGrid<value_tt,interleave>::ngr;
  using UniformGrid<value_tt,interleave>::nst;
  using UniformGrid<value_tt,interleave>::nen;
  /// memory settings
  using UniformGrid<value_tt,interleave>::mgr;
  using UniformGrid<value_tt,interleave>::mst;
  using UniformGrid<value_tt,interleave>::men;
  using UniformGrid<value_tt,interleave>::SZ2; //< mgr[1]*mgr[2], for unpack_ind
  using UniformGrid<value_tt,interleave>::GSZ;
  using UniformGrid<value_tt,interleave>::nintz; // z-size of the interleaved data
  /// iterator settings
  using UniformGrid<value_tt,interleave>::igr;
  using UniformGrid<value_tt,interleave>::ist;
  using UniformGrid<value_tt,interleave>::ien;
  
  
  
  vector<vec_type> x[3];

  mngptr<SpaceRegion> conf;

  int make_regular(){
    for(int i=0;i<3;i++){
      vec_type *g1d = new vec_type [ngr[i]];
      for(int j=0;j<ngr[i];j++)
        g1d[j]=pref[i]+j*dx[i];
      SetAxis(i,g1d,g1d+ngr[i]);
      delete[]g1d;
    }
    return 1;
  }

public:

  void init(const Vector_3 &v1, const Vector_3&v2, const int dir, const int *sz, int *dim_type_){
    UniformGrid<value_tt,interleave>::init(v1,v2,dir,sz);
    for(int i=0;i<3;i++)
      dim_type[i]=dim_type_[i];
    make_regular();
  }

  template<class inp_it>
  int SetAxis(size_t aind, inp_it beg, inp_it end){
    size_t nv=0;
    x[aind].clear();
    for(;beg!=end;++beg,++nv)
      x[aind].push_back(*beg);

    int sz[3]={ngr[0],ngr[1],ngr[2]};
    Vector_3 p1=b.get_p1(),p2=b.get_p2();
    if(nv>=2){
      sort(x[aind].begin(),x[aind].end());
      pref[aind]=x[aind][0];
      p1[aind]=x[aind][0],p2[aind]=x[aind][nv-1];
      dx[aind]=(p2[aind]-p1[aind])/(nv-1);
      sz[aind]=nv;
    }
    else{ // inactive dimension
      pref[aind] = nv>0 ? x[aind][0] : 0;
      p1[aind]=p2[aind]=pref[aind];
      dx[aind]=0;
      sz[aind]=1;
    }
    UniformGrid<value_tt,interleave>::init(p1,p2,0,sz);
    return 1;
  }

  int SetConfinement(mngarg<SpaceRegion> conf_){
    conf.reset(conf_);
    return 1;
  }

  Vector_3 Position(int ix, int iy, int iz) const{
    return Vector_3(x[0][ix],x[1][iy],x[2][iz]);
  }

  int GetCoeff(const Vector_3 &place, int *int_ind, vec_type *values, int force_external=0, vector<Vector_3> *nonlocal=NULL, bool all_nonlocal=false) const;

};

template <class value_tt>
class InterpBox:  public Box{
public:
  value_tt cube[2][2][2];
  
  /// gets the value by linear index 
  value_tt &operator[](int i) const {
    return *((value_tt *)cube+i);
  }

  int pack_ind(int i, int j, int k) const {
    return 4*i+2*j+k;
  }

  void unpack_ind(int ind, int &ix, int &iy, int &iz) const {
    ix=ind/4;
    iy=(ind-ix*4)/2;
    iz=ind%2;
  }

  InterpBox(const Vector_3 &sp1, const Vector_3 &sp2): Box(sp1,sp2){}

  /// Gets interpolation coefficients and indicies
  /// @return the number of nonzero indicies
  int GetCoeff(const Vector_3 &place, int *int_ind, vec_type *values, int force_external=0) const;
};

# if 0

struct GrdData{
 vector<Vector_3> points;
 vector<vec_type> coeffs;
 vector<int> shifts;
 int mind;
};

/// given a set of points, creates a box for finding interpolation coefficients
/// to use in 3-linear interpolation
template <class grid_it> 
int InterGridInterpolation(grid_it beg, grid_it end, const Vector_3 &place){
  vector<int> shifts;
  grid_it it=beg;
  for(;it!=end;++it){
    arr[i].mind=it->GetShifts(place,arr[i].shifts,arr[i].coeffs,1);
    // detecting the points
    int ni=arr[i].coeffs.size();
    for(j=0;j<ni;j++){
      points.push_back(it->get_position(arr[i].mind+arr[i].shifts[j]);
    }
  }
  // building the box
  Vector_3 p1, p2;
  get_min_box(points.begin(), points.end(), place,p1,p2);
  InterpBox ibox(p1,p2);
  int bind[8];
  vec_type bcoeff[8];
  int nb=ibox.GetCoeff(place,bind,bcoeff,1);
  for(i=0;i<nb;i++){
    Vector_3 grdpoint=ibox.get_position(bind[i]);
    /// detecting which grid is better to represent this particular point

  }
}

# endif
# endif

