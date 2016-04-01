/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2009        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.14 $
 *   $Date: 2014/04/17 13:51:00 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/splines.h,v 1.14 2014/04/17 13:51:00 kazeev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/splines.h,v $
$Revision: 1.14 $
$Author: kazeev $
$Date: 2014/04/17 13:51:00 $
*/
/*s****************************************************************************
 * $Log: splines.h,v $
 * Revision 1.14  2014/04/17 13:51:00  kazeev
 * Fixed gcc compilation. Added BoxHamiltonian.
 *
 * Revision 1.13  2012/09/19 16:37:25  morozov
 * Made compartible with Linux Intel compiler. Makefile is updated.
 *
 * Revision 1.12  2012/09/06 08:27:34  valuev
 * added table interaction
 *
 * Revision 1.11  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.3  2011/01/31 18:33:03  lesha
 * some comments are included
 *
 * Revision 1.2  2011/01/31 16:53:51  valuev
 * added MKL configurations
 *
 * Revision 1.1  2010/12/29 22:31:34  valuev
 * added splines
 *
 * Revision 1.8  2010/09/08 08:13:48  valuev
 * compiled with icc
 *
 * Revision 1.7  2010/06/12 17:57:31  valuev
 * some workflow coding
 *
 * Revision 1.6  2009/12/23 15:09:41  valuev
 * job manager+parametric update
 *
 * Revision 1.5  2009/07/24 05:08:46  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.4  2009/06/10 20:53:33  valuev
 * updated splindes, trajReader
 *
 * Revision 1.3  2009/06/08 08:01:05  valuev
 * updated splines
 *
 * Revision 1.2  2009/04/28 02:56:56  valuev
 * added splines derivative
 *
 * Revision 1.1  2009/04/27 18:59:40  valuev
 * Added splines
 *
*******************************************************************************/
#ifndef _SPLINES_H
#define _SPLINES_H

//e \file splines.h \brief Templates for multidimensional cubic splines.

# include <utility>
# include <vector>
# include <map>
# include <limits>
# include "refobj.h"
# include "pencil.h"
# include "logexc.h"

# ifdef USE_MKL
# include "mkl.h"
# define USE_LAPACK 1
# endif

using namespace std;

///\en Maximal dimension
#define MAX_DIM 10

//e tensor of arbitrary dimension, all elementrs stored
template <class T>
class Tensor{
protected:
  size_t dim;
  size_t sz; //e< whole data size
  mngptr<T> parr; // data array
  vector<size_t> count; // elements number along each dimension
  vector<size_t> vsz; // span for each dimension (used in indices packing)
public:
  typedef T value_t;
  typedef T value_type;

  class iterator1d {
  friend class Tensor<T>;
    size_t span;
    T* ptr;
    iterator1d(T* ptr_, size_t span_):ptr(ptr_),span(span_){}
  public:
    typedef forward_iterator_tag iterator_category;
    typedef T value_type;
    typedef int difference_type;
    typedef int distance_type;	// retained
    typedef T* pointer;
    typedef T& reference;
    //e iterator difference, compatibility is not checked
    int operator-(const iterator1d& other) const {
      int res=ptr-other.ptr;
      return res/span;
    }

    T& operator*() const {
      return *ptr;
    }

    T* operator->() const {
      return ptr;
    }


    iterator1d& operator++(){
      ptr+=span;
      return *this;
    }

    //e positive increment to an iterator
    iterator1d& operator+=(int incr){
      ptr+=span*incr;
      return *this;
    }

    iterator1d operator++(int){ // postfix
      iterator1d tmp=*this;
      ++*this;
      return tmp;
    }

    bool operator!=(const iterator1d &other) const{
      return ptr!=other.ptr;
    }

  };

  Tensor():dim(0){}

  //e constructor for specifying up to 4 dimensions, own array
  Tensor(size_t dim_, size_t cnt0, size_t cnt1=1, size_t cnt2=1,size_t cnt3=1):dim(dim_),sz(0){
    set_dim(dim);
    set_count(0,cnt0);
    set_count(1,cnt1);
    set_count(2,cnt2);
    set_count(3,cnt3);
    init(NULL,1);
  }

  //e general constructor, if arr=NULL and managed=1, allocates own array
  template <class ind_it>
  Tensor(size_t dim_, ind_it ax_counts, T *arr, int managed=1):dim(dim_),sz(0){
    set_dim(dim);
    for(size_t i=0;i<dim;i++)
      set_count(i,*ax_counts++);
    init(arr,managed);
  }
  void set_dim(size_t dim_){
    dim=dim_;
    count.resize(dim,1);
  }

  size_t dimension() const {
    return dim;
  }

  //e sets the count of values along specified axis, init must be called after all axes are setup
  void set_count(size_t ax, size_t num){
    if(ax<dim)
      count[ax]=num;
  }

  //e gets the count along given axis
  size_t get_count(size_t ax) const {
    return ax<dim ? count[ax] : 0;
  }

  int init(T *arr=NULL, int managed=1){
    vsz.resize(dim);
    sz=1;
    for(size_t i=0;i<dim;i++){
      vsz[dim-i-1]=sz;
      sz*=count[dim-i-1];
    }
    if(!arr && managed)
      parr.reset(new T[sz],1|0x8);
    else
      parr.reset(arr,managed ? managed|0x8 : 0);
    return 1;
  }

  //e sets the whole tensor
  void Set(const T &val){
    for(size_t i=0;i<sz;i++)
      parr[i]=val;
  }

  //e gets full size in memory
  size_t msize() const {
    return sz;
  }

  //e gets full size that may be addressed by packed index
  size_t size() const {
    return sz;
  }

  //e gets the data reference by packed index
  T &raw(size_t ind){
     return parr[ind];
  }


  //e gets the data reference by tensor index set
  T &operator()(size_t i0,size_t i1=0, size_t i2=0, size_t i3=0){
    size_t inda[MAX_DIM]={i0,i1,i2,i3};
    for(size_t i=4;i<MAX_DIM;i++)
      inda[i]=0;
    return get_value(inda);
  }

  //e gets the data reference by tensor index set
  template <class ind_it>
  T &get_value(ind_it inds){
    return parr[pack_ind(inds)];
  }

  //e gets the data value by tensor index set
  template <class ind_it>
  T get_value(ind_it inds) const {
    return parr[pack_ind(inds)];
  }

  template <class ind_it>
  size_t pack_base_ind(ind_it inds) const{
    size_t ind=0;
    for(size_t i=0;i<dim;i++)
      ind=ind*count[i] + *inds++;
    return ind;
  }

  template <class ind_it>
  void unpack_base_ind(size_t ind, ind_it inds) const {
    for(size_t i=0; i<dim-1;i++){
      *inds++=ind/vsz[i];
      ind%=vsz[i];
    }
    *inds=ind;
  }

  template <class ind_it>
  size_t pack_ind(ind_it inds) const{
    return pack_base_ind(inds);
  }

  template <class ind_it>
  void unpack_ind(size_t ind,ind_it inds){
    unpack_base_ind(ind,inds);
  }

  template <class ind_it>
  iterator1d axis_begin(ind_it inds,size_t axis){
    return iterator1d(parr.ptr()+pack_ind(inds),vsz[axis]);
  }

};


// tridiagonal matrix algorithm
// 3-point uniform Gaussian exclusion, a[0] and b[n-1] are ignored
// c0 b0 0  ..
// a1 c1 b1 0 ..
// 0  a2 c2 b2 0 ..
// ..
// ui is a solution input iterator
// if not NULL, working array w_arr must provide a space of n*sizeof(coeff_t)+n*sizeof(func_t)
// w_arr may intersect (keeping element alignment!) with argument iterators
template<class arg_it, class coeff_it, class func_it>
void tdma_solve(size_t n, coeff_it ai, coeff_it bi, coeff_it ci,func_it fi, arg_it ui, void *w_arr=NULL){
  int res=0;
  typedef typename iterator_traits<coeff_it>::value_type coeff_t;
  typedef typename iterator_traits<func_it>::value_type func_t;

  coeff_t *l;
  func_t *k;

  if(!w_arr){
    l= new coeff_t[n];
    k= new func_t[n];
    res=1;
  }
  else{
    l=(coeff_t *)w_arr;
    k=(func_t *)((char *)w_arr+sizeof(coeff_t)*n);
  }

  coeff_t num=*ci;
  l[0]=-(*bi)/num;
  k[0]=(*fi)/num;
  ++ci;
  ++fi;
  ++bi;
  ++ai;
  for(size_t i=1;i<n;i++,++ci,++fi,++bi,++ai){
    coeff_t num=*ci+(*ai)*l[i-1];
    if(fabs(num)<1./(1e-2*numeric_limits<coeff_t>::max())){
      l[i]=(coeff_t)0;
      k[i]=(func_t)0;
    }
    else{
      l[i]=-(*bi)/num;
      k[i]=(*fi-(*ai)*k[i-1])/num;
    }
  }
  for(int i=(int)n-2;i>=0;i--)
    k[i]=l[i]*k[i+1]+k[i];

  for(size_t i=0;i<n;i++)
    *ui++=k[i];

  if(res){
    delete [] k;
    delete [] l;
  }
}

//e flags for spline_solve and spline_matrix
enum SPLINE_FLAGS {
  SPLINE_NORMAL=0,
  SPLINE_NATURAL=1,
  SPLINE_PERIODIC=2,
  SPLINE_BOUNDARY=0xff, // mask
  SPLINE_AKIMA=0x100,
  SPLINE_TYPE=0xff00 // mask
};

//e fills the matrix with spline equations
//e depending on the filler may fill the 3-diagonal matrix and right-hand f for usual spline
//e or more complicated matix with f at the left hand as variables
//e The first and last equation are formed by special rules:
//e for NORMAL spline: first derivative at left and right end is equal to first_eq[0] and first_eq[1] respectively
//e for NATURAL spline: the first equation is c0*first_eq[0]+c1*first_eq[1]+f0*first_eq[2]+f1*first_eq[3]=first_eq[4],
//e                     the last equation is  c(n-1)*last_eq[0]+cn*last_eq[1]+f(n-1)*last_eq[2]+fn*last_eq[3]=last_eq[4]
//e note: in set_equation the free item is specified as if were standing at the right hand side
template<class arg_it, class filler_t, class cond_it>
int spline_matrix(size_t n,arg_it xi,filler_t &filler, int flag, cond_it first_eq, cond_it last_eq){
  typedef typename iterator_traits<arg_it>::value_type arg_t;
  arg_t x0=*xi++;
  arg_t x1=*xi++;
  arg_t d0, d1=x1-x0;
  //a[0]=(arg_t)0;
  if((flag&SPLINE_BOUNDARY)==SPLINE_NORMAL){
    filler.set_equation(0, 0/*a*/,1./*c*/,0./*b*/,0./*fm1*/,0./*f0*/,0./*f1*/,*first_eq/*right*/);
  }
  else if((flag&SPLINE_BOUNDARY)==SPLINE_NATURAL){
    filler.set_equation(0, 0/*a*/,*first_eq/*c*/,*(first_eq+1)/*b*/,0./*fm1*/,*(first_eq+2)/*f0*/,*(first_eq+3)/*f1*/,*(first_eq+4)/*right*/);
  }
  else{ // periodic
    return
      LOGERR(-1,"spline_matrix: periodic spline is not implemented",LINFO);
  }
  for(size_t i=1;i<n-1;i++,++xi/*,++fi*/){
    x0=x1;
    x1=*xi;
    d0=d1;
    d1=x1-x0;
    arg_t k1=d1/d0;
    if((flag&SPLINE_TYPE)==SPLINE_AKIMA){
      filler.set_equation(i, 0./*a*/,(d1+d0)/*c*/,0./*b*/,k1/*fm1*/,(-k1+1./k1)/*f0*/,-1./k1/*f1*/,0./*right*/);
    }
    else{
      //filler.set_equation(i, 0/*a*/,1/*c*/,0/*b*/,0/*fm1*/,0/*f0*/,0/*f1*/,0./*right*/);
      filler.set_equation(i, d1/*a*/,2*(d1+d0)/*c*/,d0/*b*/,3*k1/*fm1*/,3*(-k1+1./k1)/*f0*/,-3./k1/*f1*/,0./*right*/);
    }
  }
  x0=x1;
  d0=d1;
  /*fm1=f0;
  f0=f1;*/
  if((flag&SPLINE_BOUNDARY)==SPLINE_NORMAL){
    filler.set_equation(n-1, 0./*a*/,1./*c*/,0./*b*/,0./*fm1*/,0./*f0*/,0./*f1*/,*last_eq/*right*/);
  }
  else if((flag&SPLINE_BOUNDARY)==SPLINE_NATURAL){
    filler.set_equation(n-1, *last_eq/*a*/,*(last_eq+1)/*c*/,0./*b*/,*(last_eq+2)/*fm1*/,*(last_eq+3)/*f0*/,0./*f1*/,*(last_eq+4)/*right*/);
  }
  else{ // periodic
    return
      LOGERR(-1,"spline_solve: periodic spline is not implemented",LINFO);
  }
  //b[n-1]=(arg_t)0;
  return 1;
}

template<class arg_t, class func_it>
class filler_3d_t{
public:
  typedef typename iterator_traits<func_it>::value_type func_t;
protected:
  func_it fi;
  func_t f0;
  func_t f1;
  func_t fm1;
  size_t n;
public:
  func_it fi0;
  arg_t *a, *b, *c;
  func_t *ff;
  filler_3d_t(size_t n_, func_it fi_, arg_t *a_, arg_t *b_, arg_t *c_, func_t *ff_):fi(fi_),fi0(fi_),a(a_),b(b_),c(c_),ff(ff_),n(n_){}
  void set_equation(size_t i, arg_t ai, arg_t ci, arg_t bi,arg_t cfm1, arg_t cf0, arg_t cf1, func_t right){
    if(i==0){
      fi=fi0;
      fm1=0.;
      f0=*fi++;
      f1=*fi++;
    }
    else{
      fm1=f0;
      f0=f1;
      if(i<n-1)
        f1=*fi++;
      else
        f1=0.;
    }
    a[i]=ai;
    c[i]=ci;
    b[i]=bi;
    ff[i]=right-(cfm1*fm1+cf0*f0+cf1*f1);
  }

};

//e cubic spline
//e inputs:
//e n -- number of nodes, xi -- sorted argument values at nodes
//e fi -- function values at nodes
//e first_eq, last_eq -- conditions at the ends
//e flag -- spline type (enum)
//e outputs:
//e ci -- first derivatives at nodes (obtained by standard cubic spline solve requiring continuous 2nd derivative)
//e w_arr, if not NULL must provide space for 3*n*sizeof(arg_t)+n*sizeof(func_t) entries
//e The first and last equations are formed by special rules:
//e for NORMAL spline: first derivative at left and right end are equal to first_eq[0] and first_eq[1] respectively
//e for NATURAL spline: the first equation is c0*first_eq[0]+c1*first_eq[1]+f0*first_eq[2]+f1*first_eq[3]=first_eq[4],
//e                     the last equation is  c(n-1)*last_eq[0]+cn*last_eq[1]+f(n-1)*last_eq[2]+fn*last_eq[3]=last_eq[4]
//e note: in set_equation the free item is specified as if were standing at the right hand side
template<class arg_it, class func_it, class coeff_it, class cond_it>
int spline_solve(size_t n,arg_it xi,func_it fi,coeff_it ci, int flag, cond_it first_eq, cond_it last_eq, void *w_arr=NULL){
  char *rarr=NULL;
  typedef typename iterator_traits<arg_it>::value_type arg_t;
  typedef typename iterator_traits<func_it>::value_type func_t;

  arg_t *a, *b, *c;
  func_t *ff;

  if(!w_arr){
    rarr = new char[3*n*sizeof(arg_t)+n*sizeof(func_t)];
    w_arr = (void *)rarr;
  }
  a=(arg_t *)w_arr;
  b=(arg_t *)((char *)w_arr+n*sizeof(arg_t));
  c=(arg_t *)((char *)w_arr+2*n*sizeof(arg_t));
  ff=(func_t *)((char *)w_arr+3*n*sizeof(arg_t));

  filler_3d_t<arg_t, func_it> filler(n,fi,a,b,c,ff);
  // filling 3-d matrix
  spline_matrix(n,xi,filler,flag,first_eq,last_eq);

  /*
  arg_t x0=*xi++;
  arg_t x1=*xi++;
  arg_t d0, d1=x1-x0;
  func_t f0=*fi++;
  func_t f1=*fi++;
  func_t fm1;

  a[0]=(arg_t)0;
  if(flag==SPLINE_NORMAL){
    c[0]=d1/3.;
    b[0]=d1/6.;
    ff[0]=(f1-f0)/d1-val0;
  }
  else if(flag==SPLINE_NATURAL){
    c[0]=(arg_t)1.;
    b[0]=(arg_t)der0;
    ff[0]=(func_t)val0;
  }
  else{ // periodic
    return
      LOGERR(-1,"spline_solve: periodic spline is not implemented",LINFO);
  }
  for(size_t i=1;i<n-1;i++,++xi,++fi){
    x0=x1;
    x1=*xi;
    d0=d1;
    d1=x1-x0;
    fm1=f0;
    f0=f1;
    f1=*fi;
    a[i]=d0/6.;
    c[i]=(d1+d0)/3.;
    b[i]=d1/6.;
    ff[i]=(f1-f0)/d1-(f0-fm1)/d0;
  }
  x0=x1;
  d0=d1;
  fm1=f0;
  f0=f1;
  if(flag==SPLINE_NORMAL){
    a[n-1]=d0/6.;
    c[n-1]=d0/3.;
    ff[n-1]=-(f0-fm1)/d0+der1;
  }
  else if(flag==SPLINE_NATURAL){
    a[n-1]=(arg_t)der1;
    c[n-1]=(arg_t)1.;
    ff[n-1]=(func_t)val1;
  }
  else{ // periodic
    return
      LOGERR(-1,"spline_solve: periodic spline is not implemented",LINFO);
  }
  b[n-1]=(arg_t)0;*/
  // solving
  tdma_solve(n,a,b,c,ff,ci,c);
  if(rarr)
    delete [] rarr;
  return 1;
}


/// TODO: check dimensions with nnodes<3 !!!
template<class arg_t, class grid_t= Tensor<arg_t> >
class GridSpline{
  mngptr<grid_t> grid_f;  //e< values grid
  typedef typename grid_t::value_t value_t;
  typedef Tensor<value_t> cgrid_t;
  refvector<cgrid_t> grid_c;  //e< coefficients grids
  size_t ngrids; // 2^dim
  vector<size_t> axis_order; // not fixed axis go first
  Tensor<int> slices; // tensor of the dimension of fixed axis number
  vector<size_t> wslice; //e< work slice (0 for fixed axis?)
  int wmode;

  // describes mesh along some axis
  struct argdescr_t {
    arg_t xs, xe, dx;
    size_t n;
    int flag; //e< 0=regular nodes, 1= user-supplied nodes (vector<arg_t> x)
    vector<arg_t> x; // mesh coordinates for unregular mesh
    int ex[2]; //e< extrapolation types at 2 ends
    int bcd[2]; //e< degree of derivative for boundary condition at left and right ends
    value_t bcv[2]; //e< value of derivative for boundary condition at left and right ends
    int idegree; //e< interpolation degree along given axis
    argdescr_t():n(0),flag(0){
      ex[0]=ex[1]=3;
      bcd[0]=bcd[1]=3;
      bcv[0]=bcv[1]=(value_t)0;
      idegree=3;
    }
  };
  size_t dim;
  size_t spdim; // number of unfixed axis
  vector<argdescr_t> arg;
  int split; // bit flag for fixed / unfixed axis

  Tensor<int> indexer;
  vector<int> matr_grids;
  size_t msize1, msize0;
  pencil<double> pmatr;
  vector< pair<value_t,arg_t> > equations;

  friend class filler_matr_t;
  class filler_matr_t{
    GridSpline<arg_t, grid_t> *parent;
    size_t grf, grc;
    Tensor<int>::iterator1d axi, axi0;
    size_t n;
    int im1, i0, i1;
  public:
    //arg_t c[2];
    filler_matr_t(GridSpline<arg_t, grid_t> *par, size_t n_, size_t grf_, size_t grc_, Tensor<int>::iterator1d axi_, Tensor<int>::iterator1d axi0_):parent(par),grf(grf_),grc(grc_), n(n_), axi0(axi0_), axi(axi_){
      // filtering the grids
      size_t k=0;
      for(;k<parent->matr_grids.size();k++){
        if(grf==parent->matr_grids[k])
          break;
      }
      if(k>=parent->matr_grids.size()){
        LOGERR(-1,"filler_matr: unused grid is supplied as argument!!!",LINFO); // this grid is not in the list
        k=0;
      }
      grf=k;
      k=0;
      for(;k<parent->matr_grids.size();k++){
        if(grc==parent->matr_grids[k])
          break;
      }
      if(k>=parent->matr_grids.size()){
        LOGERR(-1,"filler_matr: unused grid is supplied as argument!!!",LINFO); // this grid is not in the list
        k=0;
      }
      grc=k;
    }
    void set_equation(size_t i, arg_t ai, arg_t ci, arg_t bi,arg_t cfm1, arg_t cf0, arg_t cf1, value_t left){
      if(i==0){
        //axi=axi0;
        im1=-1;
        i0=axi.operator->()-axi0.operator->();
        axi++;
        i1=axi.operator->()-axi0.operator->();
        axi++;
      }
      else{
        im1=i0;
        i0=i1;
        if(i<n-1){
          i1=axi.operator->()-axi0.operator->();
          axi++;
        }
        else
          i1=-1;
      }
      size_t eqn=parent->equations.size();
      int indf=eqn*parent->msize0+parent->indexer.size()*grf;
      int indc=eqn*parent->msize0+parent->indexer.size()*grc;
      // clearing entries
      for(size_t j=0;j<parent->msize0;j++)
        parent->pmatr[eqn*parent->msize0+j]=0;
      if(im1>=0){
        parent->pmatr[indf+im1]=cfm1;
        parent->pmatr[indc+im1]=ai;
      }
      parent->pmatr[indf+i0]=cf0;
      parent->pmatr[indc+i0]=ci;
      if(i1>=0){
        parent->pmatr[indf+i1]=cf1;
        parent->pmatr[indc+i1]=bi;
      }
      parent->equations.push_back(make_pair(left,1.));
    }
  };


  int make_spline(size_t start_lev, size_t end_lev){
    arg_t zero_eq0[]={1.,0.,0.,0.,0.};
    arg_t zero_eq1[]={0.,1.,0.,0.,0.};

    // finding maximal number of entries for w_arr
    size_t nmax=0;
    for(size_t sp=start_lev;sp<end_lev;sp++){ // spline level
      size_t curax=axis_order[sp];
      if(arg[curax].n>nmax)
        nmax=arg[curax].n;
    }
    char *w_arr= NULL;
    if(wmode==0)
      w_arr=new char[nmax*(3*sizeof(arg_t)+sizeof(value_t))];

    for(size_t sp=start_lev;sp<end_lev;sp++){ // spline level
      // prepare dimension loop counters
      vector<size_t> start, lim, counts;
      vector<arg_t> xi;
      size_t curax=axis_order[sp];
      // filling arguments
      if(arg[curax].flag)
        xi=arg[curax].x;
      else{
        arg_t xx=arg[curax].xs;
        for(size_t i=0;i<arg[curax].n;i++){
          xi.push_back(xx);
          xx+=arg[curax].dx;
        }
      }


      for(size_t i=0;i<dim;i++){
        if(i==curax){ // this is spline axis
          if(split&(1<<i))// strange, spline axis is fixed, report error
            return LOGERR(-1,fmt("GridSpline.make_spline: requested interpolation along fixed axis %d!",i),LINFO);
          start.push_back(0);
          counts.push_back(0);
          lim.push_back(1);
        }
        else if(split&(1<<i)){ // this is fixed axis
          start.push_back(wslice[i]);
          counts.push_back(wslice[i]);
          lim.push_back(wslice[i]+1);
        }
        else{ // variation axis
          lim.push_back(arg[i].n);
          start.push_back(0);
          counts.push_back(0);
        }
      }
      // dimension loop
      do{
        // checking limits
        size_t i=0;
        for(;i<dim;i++){
          if(i!=0)
            counts[i]++;
          if(counts[i]>=lim[i])
            counts[i]=start[i];
          else
            break;
        }
        if(i>=dim) // loop over
          break;

        vector<size_t> icounts, icounts0(spdim,0);
        if(wmode){ // filling the indicies for matrix adressing
          for(size_t j=0;j<spdim;j++)
            icounts.push_back(counts[axis_order[j]]);
        }

        size_t nu=1<<sp;
        for(size_t gr=0;gr<nu;gr++){ // grid loop
          size_t grc0=gr+nu, grc=0, grf=0;
          for(size_t j=0;j<dim;j++){
            if(grc0&(1<<j))
              grc|=(1<<axis_order[j]);
            if(gr&(1<<j))
              grf|=(1<<axis_order[j]);
          }

          typename cgrid_t::iterator1d c_it=grid_c[grc-1]->axis_begin(counts.begin(),curax);
          if(grf==0){ // f-grid is a special case
            // derivative degrees for not-a-node axes
            //arg_t a1,a0,b1,b0;
            arg_t first_eq[5], last_eq[5];
            // right side
            arg_t xx=arg[curax].flag ? arg[curax].x[arg[curax].n-1]-arg[curax].x[arg[curax].n-2] : arg[curax].dx;
            get_splfun(arg[curax].idegree,(size_t)arg[curax].bcd[1],0,0,xx,last_eq[2]/*a0*/,last_eq[0]/*b0*/,last_eq[3]/*a1*/,last_eq[1]/*b1*/);
            last_eq[4]=arg[curax].bcv[1];

            // left side
            xx=arg[curax].flag ? arg[curax].x[1]-arg[curax].x[0] : arg[curax].dx;
            get_splfun(arg[curax].idegree,(size_t)arg[curax].bcd[0],0,xx,xx,first_eq[2]/*a0*/,first_eq[0]/*b0*/,first_eq[3]/*a1*/,first_eq[1]/*b1*/);
            first_eq[4]=arg[curax].bcv[0];

            if(wmode==0){
              typename grid_t::iterator1d f_it=grid_f->axis_begin(counts.begin(),curax);
              spline_solve(arg[curax].n,xi.begin(),f_it,c_it,flag,first_eq,last_eq,w_arr);
            }
            else{
              Tensor<int>::iterator1d i_it=indexer.axis_begin(icounts.begin(),sp),i_it0=indexer.axis_begin(icounts0.begin(),sp);
              filler_matr_t filler(this,arg[curax].n,grf,grc,i_it,i_it0);
              spline_matrix(arg[curax].n,xi.begin(),filler,flag,first_eq,last_eq);
            }
          }
          else{
            if(wmode==0){
              //typename grid_t::iterator1d f_it=grid_f->axis_begin(counts.begin(),curax);
              typename cgrid_t::iterator1d f_it=grid_c[grf-1]->axis_begin(counts.begin(),curax);
              spline_solve(arg[curax].n,xi.begin(),f_it,c_it,SPLINE_NATURAL,zero_eq0,zero_eq1,w_arr);
            }
            else{
              Tensor<int>::iterator1d i_it=indexer.axis_begin(icounts.begin(),sp),i_it0=indexer.axis_begin(icounts0.begin(),sp);
              filler_matr_t filler(this,arg[curax].n,grf,grc,i_it,i_it0);
              spline_matrix(arg[curax].n,xi.begin(),filler,SPLINE_NATURAL,zero_eq0,zero_eq1);
            }
          }
        }
        counts[0]++;
      }while(1);

    }
    if(w_arr)
      delete [] w_arr;
    return 1;
  }

protected:
  template <class ogrid_t>
  grid_t *alloc_grid(ogrid_t *dummy, size_t dim_, size_t cnt0, size_t cnt1, size_t cnt2,size_t cnt3) const { return NULL; }

  /* Commtented by I.Morozov for compartibility with Linux Intel compiler
	template <>
  grid_t *alloc_grid(Tensor<arg_t> *dummy, size_t dim_, size_t cnt0, size_t cnt1, size_t cnt2,size_t cnt3) const {
    return new Tensor<arg_t>(dim_,cnt0,cnt1,cnt2,cnt3);
  } */

  template <class ogrid_t, class ind_it>
  grid_t *alloc_grid(ogrid_t *dummy, size_t dim_, ind_it ax_counts) const { return NULL; }

  template <class ind_it>
  grid_t *alloc_grid(Tensor<arg_t> *dummy, size_t dim_, ind_it ax_counts) const {
    return new Tensor<arg_t>(dim_,ax_counts, NULL,1);
  }


  int flag;
public:
  GridSpline():dim(0),spdim(0),split(0),ngrids(0),wmode(0),grid_c(1),flag(SPLINE_NATURAL){}

  //e own array, simplified for grid_t=Tensor, no SetFuncGrid needed
  GridSpline(size_t dim_, size_t cnt0, size_t cnt1=1, size_t cnt2=1,size_t cnt3=1):dim(0),spdim(0),split(0),ngrids(0),wmode(0),grid_c(1),flag(SPLINE_NATURAL){
    SetFuncGrid(alloc_grid(grid_f.ptr(),dim_,cnt0,cnt1,cnt2,cnt3),1);
  }


  void SetFuncGrid(grid_t *grdf=NULL, int mangrdf=0){
    dim=0; // indicates the new construction
    grid_f.reset(grdf,mangrdf);
  }

  //e Set the spline  and boundary condition type (see SPLINE_FLAGS enum)
  int SetMode(int flag_){
    swap(flag,flag_);
    return flag_;
  }
  // split bit flag, bit==1 means equation direction
  int BeginBuild(int splitflag=0){
    dim=grid_f->dimension();
    wslice.resize(dim);
    arg.resize(dim);
    for(size_t i=0;i<dim;i++){ //default settings: all dimensions have range 0:1
      size_t nv =grid_f->get_count(i);
      arg[i].n=nv;
      arg[i].xs=(arg_t)0.;
      arg[i].xe=nv >=2 ? (arg_t)1. : arg[i].xs;
      arg[i].dx=nv >=2 ? (arg_t)(1./(nv-1)) : (arg_t) 0;
      arg[i].flag=0;

    }
    ngrids=1<<dim;
    grid_c.clear();
    //grid_c.resize(ngrids-1);
    for(size_t i=0;i<ngrids-1;i++){
      Tensor<value_t> *grd=new Tensor<value_t>();
      grid_c.push_back(grd);
      grd->set_dim(dim);
      for(size_t j=0; j<dim; j++)
        grd->set_count(j,grid_f->get_count(j));
      if(grd->init()<0)
        return LOGERR(-1,fmt("GridSpline.BeginBuild: coefficients grid[%d] initialization failed.",i),LINFO);
      grd->Set((arg_t)0);
    }
    split=splitflag;
    // forming the axis order: fixed at the end
    axis_order.resize(dim);
    size_t ia=0;
    spdim=0;
    for(size_t i=0;i<dim;i++){
      if(!(split&(1<<i))){ // this is NOT the fixed axis
        axis_order[ia++]=i;
        spdim++;
      }
    }

    size_t nfixed=dim-spdim;
    if(nfixed>0)
      slices.set_dim(nfixed);

    size_t ax=0;
    for(size_t i=0;i<dim;i++){
      if(split&(1<<i)){ // this is the fixed axis
        axis_order[ia++]=i;
        slices.set_count(ax++,arg[i].n);
      }
      wslice[i]=0;
    }
    if(nfixed>0){
      slices.init();
      slices.Set(0);
    }
    return 1;
  }

  int SetAxis(size_t aind, arg_t xs=0., arg_t xe=1.){
    if(aind>=dim)
      return LOGERR(-1,fmt("GridSpline.SetAxis: invlaid axis number (%d)",aind),LINFO);
    size_t nv=arg[aind].n;
    arg[aind].xs=xs;
    arg[aind].xe=xe;
    arg[aind].x.clear();
    arg[aind].dx= nv>=2 ? (xe-xs)/(nv-1): (arg_t)0.;
    arg[aind].flag=0;
    return 1;
  }

  int CommitGrdAxis(size_t aind){
    if(arg[aind].flag){
      sort(arg[aind].x.begin(),arg[aind].x.end());
      arg[aind].xs=arg[aind].x[0];
      arg[aind].xe=arg[aind].x[arg[aind].n-1];
      arg[aind].dx= (arg[aind].xe- arg[aind].xs)/(arg[aind].n-1);
    }
    return 1;
  }

  template<class inp_it>
  int SetAxis(size_t aind, inp_it beg, inp_it end){
    if(aind>=dim)
      return LOGERR(-1,fmt("GridSpline.SetAxis: invlaid axis number (%d)",aind),LINFO);
    size_t nv=0;
    arg[aind].x.clear();
    for(;beg!=end;++beg,++nv)
      arg[aind].x.push_back(*beg);

    if(nv<arg[aind].n)
      return LOGERR(-1,fmt("GridSpline.SetAxis: the number of given nodes (%d) for axis %d is less than required for the grid (%d)",nv,aind,arg[aind].n),LINFO);

    if(nv>=2){
      sort(arg[aind].x.begin(),arg[aind].x.end());
      arg[aind].xs=arg[aind].x[0];
      arg[aind].xe=arg[aind].x[nv-1];
      arg[aind].dx=(arg[aind].xe-arg[aind].xs)/(nv-1);
      arg[aind].flag=1;
    }
    else{ // inactive dimension
      arg[aind].xs= nv>0 ? arg[aind].x[0] :(arg_t) 0.;
      arg[aind].xe=arg[aind].xs;
      arg[aind].dx= (arg_t)0.;
      arg[aind].flag=0;
    }
    return 1;
  }

  arg_t GetAxisMin(size_t aind) const {
    return aind<dim ? arg[aind].xs : 0;
  }

  arg_t GetAxisMax(size_t aind) const {
    return aind<dim ? arg[aind].xe : 0;
  }

  //e mode=0: direct (TDMA), =1: equation
  //e iterator slice sets the indicies for fixed axes
  template< class ind_it >
  int BeginSubgrid(ind_it slice, int mode=0){
    if(split)
      slices.get_value(slice)=1;
    wmode=mode;
    for(size_t i=0;i<dim;i++){
      if(split&(1<<i)){ // this is the fixed axis
        wslice[i]=*slice++;
      }
      else
        wslice[i]=0; // this is working axis
    }

    if(wmode){ //matrix
      // calculating the list of matrix grids and matrix dimension
      matr_grids.clear();
      size_t nu=1<<spdim;
      for(size_t gr=0;gr<nu;gr++){ // grid loop
        size_t grc=0;
        for(size_t j=0;j<dim;j++){
          if(gr&(1<<j))
            grc|=(1<<axis_order[j]);
        }
        matr_grids.push_back(grc);
      }

      indexer.set_dim(spdim);
      size_t nax=1;
      for(size_t j=0;j<spdim;j++){
        indexer.set_count(j,grid_f->get_count(axis_order[j]));
      }
      indexer.init(NULL,0); // fake tensor just to get indicies
      // allocating matrix
      msize0=msize1=matr_grids.size()*indexer.size();
      pmatr.resize(msize0*msize1,msize0*msize1);
      equations.clear();
      return make_spline(0,spdim); // this will fill default equations
    }
    return 1;
  }

  // WARNING(kazeevn) gcc can't find ant msize
  /* size_t GetMatrixSize() const { */
  /*   return wmode ? msize() : 0; */
  /* } */

  template< class ind_it >
  value_t &grd_val(ind_it ind){
    //refine(ind);
    return grid_f->get_value(ind);
  }

  value_t &grd_value(size_t i0,size_t i1=0, size_t i2=0, size_t i3=0){
    size_t inda[4]={i0,i1,i2,i3};
    return grd_val(inda);
  }

  //e set the equation c0*x0+c1*x1+...=val
  //e xi contains dim entries: x0={*(xi+0),*(xi+1),....,*(xi+dim-1)},
  //e                          x1={*(xi+dim),*(xi+dim+1),....,*(xi+2*dim-1)},...
  //e values for the fixed axes are taken from the setup slice
  template< class arg_it, class coeff_it>
  int AddEquation(size_t nterms, arg_it xi, coeff_it ci, value_t val, arg_t weight=1.){
    if(!nterms)
      return -1;
    if(!wmode)
      return -2;
    vector<arg_t> xx(dim), xx0(dim), xx1(dim);
    vector<size_t> nn(dim,0);
    vector<size_t> ind(dim);
    size_t eqn=equations.size();
    if(eqn>=msize1){ // preparing least squares matrix
      msize1=2*msize1;
      pmatr.resize(msize0*msize1,msize0*msize1);
    }
    for(size_t i=0;i<msize0;i++) // clearing matrix array
      pmatr[eqn*msize0+i]=0.;
    for(size_t i=0;i<nterms;i++){
      // forming the arg
      for(size_t j=0;j<dim;j++){
        if(j<spdim)
          xx[axis_order[j]]=*xi;
        else
          xx[axis_order[j]]=arg[axis_order[j]].x[wslice[j]];
        ++xi;
      }
      get_ind(xx.begin(),ind.begin(),xx0.begin(),xx1.begin());
      load_coeff_t lc(this);
      cell_eval(nn.begin(),xx.begin(),ind,xx0,xx1,lc);
      // packing the matrix
      for(size_t k=0;k<lc.ind.size();k++)
        pmatr[eqn*msize0+indexer.size()*lc.grd[k]+lc.ind[k]]+=lc.c[k]*(*ci);
      ++ci;
    }
    equations.push_back(make_pair(val,weight));
    return (int)equations.size()-1;
  }

  int EndSubgrid(){
    if(wmode==0) //makes spline for the selected subgrid
      return make_spline(0,spdim);
    else
      return solve_matrix();
  }

  int solve_matrix(){
    pencil<double> sqm;
    pencil<double> righth(new arg_t[msize0],msize0,msize0,1);

    // analyzing the number of equations
    size_t eqn=equations.size();
    if(eqn<msize0)
      return -(int)(msize0-equations.size()); // need more equations
    else if(eqn>msize0){ // overcomplete system, making least squares
      sqm.resize(msize0*msize0,msize0*msize0);
      for(size_t i=0;i<msize0;i++){
        for(size_t j=0;j<msize0;j++){
          sqm[msize0*i+j]=0.;
          for(size_t k=0;k<eqn;k++)
            sqm[msize0*i+j]+=pmatr[msize0*k+i]*pmatr[msize0*k+j]*equations[k].second;
        }
        righth[i]=0;
        for(size_t k=0;k<eqn;k++)
          righth[i]+=pmatr[msize0*k+i]*equations[k].first*equations[k].second;

      }
    }
    else{ // trying to solve one-to one equation system
      sqm=pmatr;
      for(size_t k=0;k<msize0;k++)
        righth[k]=equations[k].first;
    }
    //solving
    /*
    FILE *f=fopen("matr.txt","wt");
    for(size_t i=0;i<msize0;i++){
      for(size_t j=0;j<msize0;j++){
        fprintf(f,"%15g ",sqm[i*msize0+j]);
      }
      fprintf(f,"= %15g\n",righth[i]);
    }
    fclose(f);*/

# ifdef USE_LAPACK
    int nrhs=1;
    int info=0;
    pencil<int> ipiv(new int[msize0],msize0,msize0,1);
    int msz=(int)msize0;
    DGETRF(&msz,&msz,sqm.get_ptr(),&msz,ipiv.get_ptr(),&info);
    if(info<0)  // can't solve, need more equations
      return -1;
    DGETRS("T",&msz,&nrhs,sqm.get_ptr(),&msz,ipiv.get_ptr(),righth.get_ptr(),&msz,&info);
    // getting the solution
    vector<size_t> uind(spdim),rind(dim);
    for(size_t i=0;i<msize0;i++){
      size_t grd=matr_grids[i/indexer.size()];
      size_t pind=i%indexer.size();
      indexer.unpack_ind(pind,uind.begin());
      size_t j=0;
      for(;j<spdim;j++)
        rind[axis_order[j]]=uind[j];
      for(;j<dim;j++)
        rind[axis_order[j]]=wslice[axis_order[j]];
      if(!grd){
        grid_f->get_value(rind.begin())=righth[i];
        //printf("f[%d]=%g\n",rind[0],righth[i]);
      }
      else{
        grid_c[grd-1]->get_value(rind.begin())=righth[i];
        //printf("c[%d]=%g\n",rind[0],righth[i]);
      }
    }
# else // insert linsys solution here?
    return LOGERR(-1,"GridSpline.solve_matrix: need a matrix inversion routine to proceed, switch USE_LAPACK compiler option on!",1);
# endif
    return 1;
  }

  int EndBuild(int clear_matrix=1){ // finishes the levels of the fixed axes
    wmode=0;
    if(clear_matrix){
      equations.clear();
      pmatr.resize(0,0);
    }
    if(split){ // completing the fixed axes
      // completeing unset slices
      size_t sz=slices.msize();
      vector<size_t> inds( slices.dimension());
      for(size_t i=0;i<sz;i++){
        if(!slices.raw(i)){
          slices.unpack_base_ind(i,inds.begin());
          BeginSubgrid(inds.begin(),0);
          EndSubgrid();
        }
      }
      split=0; // now the unfixed axes go
      return make_spline(spdim,dim);
    }
    else
      return make_spline(0,dim); // making the whole grid
  }

  arg_t get_arg(size_t ax, size_t num) const {
    if(arg[ax].flag)
      return arg[ax].x[num];
    else{
      return arg[ax].xs+num*arg[ax].dx;
    }
  }

  // if axis is regular it becomes irregular
  bool set_arg(size_t ax, size_t num, const arg_t &val){
    if(ax>=dim)
      return false;
    if(num>=arg[ax].n)
      return false;
    if(!arg[ax].flag){ // resetting to array
      arg[ax].x.clear();
      for(size_t i=0;i<arg[ax].n;i++)
        arg[ax].x.push_back(arg[ax].xs+i*arg[ax].dx);
      arg[ax].flag=1;
    }
    arg[ax].x[num]=val;

    return true;
  }

  //e set extrapolation to the left and right of the arguments range along given axis
  //e pol_left, pol_right are polynomial degrees of the extrapolating function:
  //e <0: zero, 0: constant, 1: linear, ... 3: cubic (the same as the boundary spline segment)
  int SetExtrapolation(size_t ax,int pol_left=3, int pol_right=3){
    if(ax>=dim)
      return -1;
    arg[ax].ex[0]=pol_left;
    arg[ax].ex[1]=pol_right;
    return 1;
  }

  //e Set interpolation degree along given axis
  //e 0= step function (step in half interval)
  //e 1= linear
  //e 2= not supported
  //e 3= cubic
  int SetInterpolation(size_t ax, int degree=3){
    if(ax>=dim)
      return -1;
    arg[ax].idegree=degree;
    return 1;
  }


  //e set not-a-node boundary conditions at the left and right end of the given axis:
  //e the derivative of the given degree (1-3) is equal to given value
  int SetEnds(size_t ax,int deg_left=3,int deg_right=3, const value_t &val_left=(value_t)0, const value_t &val_right=(value_t)0){
    if(ax>=dim)
      return -1;
    if(deg_left<1 || deg_left>3)
      return -2;
    if(deg_right<1 || deg_right>3)
      return -2;

    arg[ax].bcd[0]=deg_left;
    arg[ax].bcd[1]=deg_right;
    arg[ax].bcv[0]=val_left;
    arg[ax].bcv[1]=val_right;
    return 1;
  }


  template< class arg_it, class ind_iit, class pnt_iit>
  void get_ind(arg_it x, ind_iit ind, pnt_iit x0, pnt_iit x1) const {
    for(size_t i=0;i<dim;i++,++x){
      if(*x<=arg[i].xs){
        *x0++=arg[i].xs;
        if(arg[i].flag)
          *x1++=arg[i].x[1];
        else
          *x1++=arg[i].xs+arg[i].dx;
        //*ind++=1;
        *ind++=0; // flag to use extrapolation
        continue;
      }
      else if(*x>=arg[i].xe){
        if(arg[i].flag)
          *x0++=arg[i].x[arg[i].n-2];
        else
          *x0++=arg[i].xe-arg[i].dx;
        *x1++=arg[i].xe;
        //*ind++=arg[i].n-1;
        *ind++=arg[i].n; // flag to use extrapolation
        continue;
      }
      if(arg[i].flag){ // all points
        typename vector<arg_t>::const_iterator lb=lower_bound(arg[i].x.begin(),arg[i].x.end(),*x);
        int ib=lb-arg[i].x.begin();
        if(ib>=(int)arg[i].n-1)
          ib=arg[i].n-1;
        else if(ib<=0)
          ib=1;
        *ind++=(size_t)ib;
        *x0++=arg[i].x[ib-1];
        *x1++=arg[i].x[ib];
      }
      else{ // just a range
        int ii=(int)((*x-arg[i].xs)/arg[i].dx);
        *ind++=ii+1;
        arg_t xx=arg[i].xs+ii*arg[i].dx;
        *x0++=xx;
        *x1++=xx+arg[i].dx;
      }
    }
  }


protected:
  friend struct sum_t;
  //e summation of coefficients into the value
  struct sum_t{
    GridSpline<arg_t, grid_t> const *parent;
    value_t sum;
    sum_t(GridSpline<arg_t, grid_t> const * par):parent(par),sum((value_t)0){};
    template<class ind_it>
    void eval(size_t gr, ind_it ind, arg_t coeff){
      if(coeff==0.)
        return;
      value_t v;
      if(gr==0)
        v=parent->grid_f->get_value(ind);
      else
        v=parent->grid_c[gr-1]->get_value(ind);
      sum+=v*coeff;
    }
  };

  friend struct load_coeff_t;
  //e loading coefficients
  struct load_coeff_t{
    GridSpline<arg_t, grid_t> const *parent;
    vector<int> ind;
    vector<size_t> grd;
    vector<arg_t> c;
    //arg_t c[2];
    load_coeff_t(GridSpline<arg_t, grid_t> const * par):parent(par){}
    template<class ind_it>
    void eval(size_t gr, ind_it uind, arg_t coeff){
      // filtering the grids
      size_t k=0;
      for(;k<parent->matr_grids.size();k++){
        if(gr==parent->matr_grids[k])
          break;
      }
      if(k>=parent->matr_grids.size())
        return; // this grid is not in the list
      // filtering indices
      vector<size_t> iind(parent->spdim);
      for(size_t i=0;i<parent->spdim;i++)
        iind[i]=uind[parent->axis_order[i]];
      ind.push_back(parent->indexer.pack_ind(iind.begin()));
      grd.push_back(k);
      c.push_back(coeff);
    }
  };

  void get_splfun(int idegree, size_t n,arg_t x0, arg_t x, arg_t x1, arg_t &a0, arg_t &b0, arg_t &a1, arg_t &b1) const {
    arg_t d=x1-x0;
    arg_t l1=x1-x;
    arg_t d2=d*d;
    arg_t l0=x-x0;
    if(idegree>=2){
      //arg_t d3=d2*d;
      if(n==0){
        arg_t l00=l0*l0/d2;
        arg_t l11=l1*l1/d2;
        b1=l0*l11;
        b0=-l00*l1;
        a1=2*b1/d+l11; //a1=l1*l1*(2*l0+d)/d3;
        a0=-2*b0/d+l00;  // a0=l0*l0*(2*l1+d)/d3;
      }
      else if(n==1){
        arg_t l10=2*l1*l0/d2;
        b1=l1*l1/d2-l10;
        b0=l0*l0/d2-l10;
        a1=2*(b1-l1/d)/d;
        a0=2*(-b0+l0/d)/d;
      }
      else if(n==2){
        b1=2*(l0-2*l1)/d2;
        b0=2*(2*l0-l1)/d2;
        a1=2*((2*l0-4*l1)/d+1.)/d2;
        a0=2*((2*l1-4*l0)/d-1.)/d2;
      }
      else if(n==3){
        b1=b0=6/d2;
        a1=2*b1/d;
        a0=-a1;
      }
      else
        a1=a0=b1=b0=(arg_t)0;
    }
    else if(idegree==1){ // linear
      if(n==0){
        b1=b0=0.;
        a0=l0/d;
        a1=l1/d;
      }
      else if(n==1){
        b1=b0=0.;
        a0=1./d;
        a1=-a0;
      }
      else {
        a0=a1=b1=b0=0.;
      }
    }
    else { // steps
      if(n==0){
        b1=b0=0.;
        if(l0>d/2){
          a0=0.;
          a1=1.;
        }
        else{
          a1=0.;
          a0=1.;
        }
      }
      else {
        a0=a1=b1=b0=0.;
      }
    }
  }

  //e flag=0: left, flag=1 right
  bool extrapol(int idegree, size_t n, int flag,int type,arg_t x, arg_t x0, arg_t x1, arg_t &a0, arg_t &a1, arg_t &b0, arg_t &b1) const {
    a0=b0=a1=b1=(arg_t)0;
    if(type<0) // zero continuation
      return true;
    else if(type<3){ // polynomial of degree less than 3
      size_t ne=(size_t)type;
      arg_t r=1;

      for(size_t ni=0;ni<=ne;ni++){
        if(n>ni)
          continue;
        size_t s=1;
        for(size_t i=ni;i!=ni-n;i--)
          s*=i;

        arg_t aa0, aa1, bb0, bb1;
        get_splfun(idegree,ni,x0,flag ? x1 : x0,x1,aa0,bb0,aa1,bb1);
        size_t k = ni>1 ? ni+1 : 1;
        a0+=aa0*s*r/k;
        b0+=bb0*s*r/k;
        a1+=aa1*s*r/k;
        b1+=bb1*s*r/k;
        r*= flag? (x-x1) : (x-x0);
      }
      return true;
    }
    return false; // means that normal spline continuation should be used (degree 3)
  }

  template<class der_it, class arg_it, class eval_t>
  int cell_eval(der_it n, arg_it x, vector<size_t> &ind, vector<arg_t> &x0, vector<arg_t> &x1, eval_t &ev) const {
    vector<size_t> inda(dim), indn(dim);
    vector<arg_t> a1(dim), a0(dim), b1(dim), b0(dim);
    for(size_t i=0;i<dim;i++,++x,++n){
      if(ind[i]<1){ // extrapolate left end
        indn[i]=1;
        if(extrapol(arg[i].idegree,*n,0,arg[i].ex[0],*x,x0[i],x1[i],a0[i],a1[i],b0[i],b1[i]))
          continue;
      }
      else if(ind[i]>=arg[i].n){ // extrapolate right end
        indn[i]=arg[i].n-1;
        if(extrapol(arg[i].idegree,*n,1,arg[i].ex[1],*x,x0[i],x1[i],a0[i],a1[i],b0[i],b1[i]))
          continue;
      }
      else
        indn[i]=ind[i];

      get_splfun(arg[i].idegree,*n,x0[i],*x,x1[i],a0[i],b0[i],a1[i],b1[i]);
    }

    size_t nu=1<<dim;
    value_t val=0.;
    for(size_t gr=0;gr<nu;gr++){
      for(size_t s=0;s<nu;s++){
        arg_t r=1.;
        for(size_t i=0;i<dim;i++){
          if(s&(1<<i)){ // right side of grid, left side of coeff
            inda[i]=indn[i];
            if(gr&(1<<i)) // c-frid
              r*=b0[i];
            else // f-grid
              r*=a0[i];
          }
          else{// left side of grid, right side of coeff
            inda[i]=indn[i]-1;
            if(gr&(1<<i)) // c-frid
              r*=b1[i];
            else // f-grid
              r*=a1[i];
          }
        }
        ev.eval(gr,inda.begin(),r);
      }
    }
    return 1;
  }
public:
  //e save the spline into ASCII file
  int Save(const char *filename, const char *form_arg="%g ", const char *form_val="%g ") const {
    FILE *f=fopen(filename,"wt");
    if(!f)
      return LOGERR(-1,fmt("GridSpline.Save: can't open file '%s' for writing!",filename),LINFO);
    fprintf(f,"# dim: %d\n",(int)dim);
    for(size_t i=0;i<dim;i++){
      fprintf(f,"# axis%d: %d %d %d %d %d\n",(int)i,(int)arg[i].n,arg[i].flag,arg[i].ex[0],arg[i].ex[1], arg[i].idegree);
    }
    fprintf(f,"#");
    for(size_t i=0;i<dim;i++)
      fprintf(f,"%d-arg[%d] ",(int)i,(int)i);
    size_t nu=1<<dim;
    for(size_t gr=0;gr<nu;gr++)
      fprintf(f,"%d-f[%X] ",(int)(dim+gr),gr);
    fprintf(f,"\n");
    vector<size_t> count(dim,0);
    do{
      for(size_t i=0;i<dim;i++){
        arg_t x=arg[i].flag ? arg[i].x[count[i]] : arg[i].xs+count[i]*arg[i].dx;
        fprintf(f,form_arg,x);
      }
      for(size_t gr=0;gr<nu;gr++){
        value_t v= gr==0 ? grid_f->get_value(count.begin()) : grid_c[gr-1]->get_value(count.begin());
        fprintf(f,form_val,v);
      }
      fprintf(f,"\n");
      for(size_t i=0;i<dim;i++){
        count[i]++;
        if(count[i]>=arg[i].n){
          if(i<dim-1)
            count[i]=0;
          if(i==0)
            fprintf(f,"\n");
        }
        else
          break;
      }
    }while(count[dim-1]<arg[dim-1].n);
    fclose(f);
    return 1;
  }

  //e Load the spline from ASCII file.
  //e The existing spline data is replaced by that from the file.
  int Load(const char *filename, const char *form_arg="%lf ", const char *form_val="%lf "){
    FILE *f=fopen(filename,"rt");
    if(!f)
      return LOGERR(-1,fmt("GridSpline.Load: can't open file '%s' for reading!",filename),LINFO);

    try{
      int lnum=1;
      char form[200];
      if(fscanf(f,"# dim: %u\n",&dim)!=1)
        throw lnum;
      vector<size_t> count(dim);
      vector<int> flag(dim), ex0(dim), ex1(dim), idegree(dim);
      char buff[2000];
      for(size_t i=0;i<dim;i++){
        fgets(buff,2000,f);
        sprintf(form,"# axis%d: %%u %%d %%d %%d %%d",(int)i);
        lnum++;
        ex0[i]=ex1[i]=idegree[i]=3;
        if(sscanf(buff,form,&count[i],&flag[i],&ex0[i],&ex1[i],&idegree[i])<2)
          throw lnum;

      }
      // constructing the grid
      SetFuncGrid(alloc_grid(grid_f.ptr(),dim,count.begin()),1);
      BeginBuild();

      lnum++;
      while(fgetc(f)!='\n')
        if(feof(f))
          throw lnum;

      size_t nu=1<<dim;
      vector<arg_t> xx(dim), xs(dim,numeric_limits<arg_t>::max()), xe(dim,numeric_limits<arg_t>::min());
      count.assign(dim,0);
      vector<vector<arg_t> > xp(dim);
      for(size_t i=0;i<dim;i++){
        if(flag[i])
          xp[i].resize(arg[i].n);
      }

      do{
        lnum++;
        for(size_t i=0;i<dim;i++){
          if(fscanf(f,form_arg,&xx[i])!=1)
            throw lnum;
        }
        for(size_t gr=0;gr<nu;gr++){
          value_t *v= gr==0 ? &(grid_f->get_value(count.begin())) : &(grid_c[gr-1]->get_value(count.begin()));
          if(fscanf(f,form_val,v)!=1)
            throw lnum;
        }
        for(size_t i=0;i<dim;i++){
          bool fill=true;
          for(size_t j=i+1;j<dim;j++){
            if(count[j]!=0){
              fill=false;
              break;
            }
          }
          if(fill){
            if(flag[i]) // fill the argument
              xp[i][count[i]]=xx[i];
            // find min-max
            if(xs[i]>xx[i])
              xs[i]=xx[i];
            if(xe[i]<xx[i])
              xe[i]=xx[i];
          }

          count[i]++;
          if(count[i]>=arg[i].n){
            if(i<dim-1)
              count[i]=0;
          }
          else
            break;
        }
      }while(count[dim-1]<arg[dim-1].n);
      fclose(f);
      // sorting the arguments
      for(size_t i=0;i<dim;i++){
        if(!flag[i])
          SetAxis(i,xs[i],xe[i]);
        else
          SetAxis(i,xp[i].begin(),xp[i].end());
        arg[i].ex[0]=ex0[i];
        arg[i].ex[1]=ex1[i];
        arg[i].idegree=idegree[i];
      }
    }
    catch(int res){
      fclose(f);
      return LOGERR(-1,fmt("GridSpline.Load: error reading spline data at line %d, file '%s'!",res,filename),0);
    }
    return 1;
  }


  template< class arg_it>
  value_t operator()(arg_it x) const {
    size_t inda[MAX_DIM];
    for(size_t i=0;i<MAX_DIM;i++)
      inda[i]=0;
    return (*this)(inda,x);
  }
  //e n-th order derivative in each direction
  template<class der_it, class arg_it>
  value_t operator()(der_it n, arg_it x) const {
    vector<size_t> ind(dim);
    vector<arg_t> x1(dim), x0(dim);
    get_ind(x,ind.begin(),x0.begin(),x1.begin());
    // here interprete ind
    sum_t s(this);
    cell_eval(n,x,ind,x0,x1,s);
    return s.sum;
  }

  value_t value(arg_t x0,arg_t x1=0.,arg_t x2=0.,arg_t x3=0.) const {
    arg_t xx[4]={x0,x1,x2,x3};
    return (*this)(xx);
  }

  value_t der(size_t n0, arg_t x0, size_t n1=0, arg_t x1=0., size_t n2=0, arg_t x2=0., size_t n3=0, arg_t x3=0.) const {
    arg_t xx[4]={x0,x1,x2,x3};
    size_t nn[4]={n0,n1,n2,n3};
    return (*this)(nn,xx);
  }


};



#endif
