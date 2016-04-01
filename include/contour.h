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
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/contour.h,v 1.6 2012/06/29 10:50:12 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/contour.h,v $
$Revision: 1.6 $
$Author: valuev $
$Date: 2012/06/29 10:50:12 $
*/
/*s****************************************************************************
 * $Log: contour.h,v $
 * Revision 1.6  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.48  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.47  2012/03/21 17:01:53  lesha
 * documentation
 *
 * Revision 1.46  2012/02/17 00:10:53  lesha
 * PtrContour is added
 *
*******************************************************************************/
#ifndef CONTOUR_H
#define CONTOUR_H

/// \en @file contour.h \brief Some useful interators (iterator at points, edges, planes)
/// and classes for 2D and 3D contours.
/// There are different contour classes that can be more convinient in different situations

#include "vector_3.h"
#include "plane_3.h"
#include "utiltl.h"
#include <vector>

using namespace std;

// records opposite coordinates of maximal rectangular parallelepiped which contains selected points
// returns center of this parallelepiped
template<class point_it,int N>
Vector_Nt<vec_type,N> GetBoundingBox(point_it it, point_it end,Vector_Nt<vec_type,N>* cube1,Vector_Nt<vec_type,N>* cube2);

///\en class to store the edge with some useful functions
///\ru возращает полезные функции ребра по двум образующим его вершинам
template <class point_it>
class edge_t{
  const point_it p1, p2; // two vertices of the edge
public:

  typedef typename point_it::value_type vec_t;

  // sp1, sp2 are vertices of the edge
  edge_t(const point_it &sp1, const point_it &sp2):p1(sp1),p2(sp2){}
  const vec_t &get_p1() const{return *p1;}
  const vec_t &get_p2() const{return *p2;}
  
  vec_t rel() const{
    return *p2-*p1;
  }

  vec_t dir(vec_type *len=NULL) const{
    vec_t v=rel();
    vec_type l=v.normalize();
    if(len)*len=l;
    return v;
  }

  vec_type len() const{
    return (p2-p1).norm();
  }
};

///\en iterator at edges of the contour, defined by iterator at the start and end vertex
///\ru итератор по ребрам контура, построенный на основе его начального и конечного итератора по его вершинам
template <class  point_it>
class generic_edge_it{
  point_it beg, end; // start and end iterator at the contour vertex
  point_it cur, cur1; // iterator at two adjoint vertices forming an edge

  void advance_cur1(){
    cur1++;
    if(!(cur1!=end))cur1=beg; // looping the second pointer
  }

public:

  typedef typename point_it::value_type vec_t;
  typedef edge_t<point_it> edge;

  // sbeg, send are start and end iterator at the contour vertex,
  // endpoint is 1 if edge iterator is end iterator
  generic_edge_it(point_it sbeg, point_it send, int endpoint=0): beg(sbeg), end(send), cur(sbeg),cur1(sbeg){
    if(!endpoint)advance_cur1();
  }

  generic_edge_it(const generic_edge_it &other): end(other.end), beg(other.beg), cur(other.cur),cur1(other.cur1){}
 
  edge operator*() const {
    return edge(cur,cur1);
  }
  
  generic_edge_it& operator++(){ // prefix
    cur++;
    advance_cur1();
    return *this;
  }

  generic_edge_it operator++(int){ // postfix
    generic_edge_it tmp = *this;
    ++*this;
    return tmp;
  }

  bool operator!=(const generic_edge_it &other) const{
    return cur!=other.cur;
  }
};

///\en iterator at planes normal to the contour plane and passing throw edges of the contour
/// Used in ProjectSimplex only.
///\ru итератор по плоскостям, проходящим через ребра контура и нормальным к плоскости контура. Используется только в ProjectSimplex.
template<class edge_it>
class normplanes_it{
  edge_it cur; // current edge
  Vector_3 n; // normal to the contour
public:

  normplanes_it(const Vector_3 &sn, edge_it beg):cur(beg),n(sn){
    n.normalize();
  }
  normplanes_it(const normplanes_it &other):cur(other.cur),n(other.n){}

  Plane_3 operator*() const {
    typename edge_it::edge cedge=*cur;
    return Plane_3(n%cedge.dir(),cedge.get_p1());
  }
  
  normplanes_it& operator++(){ // prefix
    cur++;
    return *this;
  }

  normplanes_it operator++(int){ // postfix
    normplanes_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const normplanes_it &other) const{
    return cur!=other.cur;
  }
};

/// Basic contour interface
template<class point_it, class edge_it=generic_edge_it<point_it>, 
  class edge_t=typename generic_edge_it<point_it>::edge>
class Contour{
public:

  typedef typename point_it::value_type vector_t;
  typedef point_it point_iterator;
  static const int dimension = vector_t::dimension;
  typedef edge_it edge_iterator;
  typedef edge_t  edge_type;

  // убрать виртуальность?
  virtual point_it points_begin() const =0;
  virtual point_it points_end() const =0;

  edge_iterator edges_begin() const{
    return edge_it(points_begin(),points_end());
  }
  edge_iterator edges_end() const{
    return edge_it(points_end(),points_end(),1); // may be (points_end(),points_begin()) ?
  }

  int GetNPoints() const {
    return 0;
  }

  int GetNEdges()const {
    return GetNPoints()-1;
  }

  /// gets vertex with given index or infinite vector if index not found
  vector_t GetPoint(int ind) const {
    point_it it=points_begin(), e=points_end();
    for(int i=0; it!=e; ++it, i++)
      if(i==ind)return *it;
    return vector_t(VEC_INFTY);
  }

  /// sets the value for given vertex
  /// returns <0 if this is not possible
  int SetPoint(int ind, const vector_t &vec){
    return -1;
  }

  void Clear(){}

  /// calculates the area of the contour (assumed flat)
  vec_type Area() const;

  /// gets contour center
  vector_t GetCenter() const;

  vector_t GetBoundingBox(vector_t *v1, vector_t *v2) const{
    return ::GetBoundingBox(points_begin(),points_end(),v1,v2);
  }
};

template<class point_it, class vector_t=typename point_it::value_type, class edge_it=generic_edge_it<point_it>, 
  class edge_t=typename generic_edge_it<point_it>::edge>
class Contour_N: public Contour<point_it,edge_it,edge_t>{};

// Contour_N specialization for 2D case
template<class point_it, class edge_it, class edge_t>
class Contour_N<point_it,Vector_2,edge_it,edge_t>: public Contour<point_it,edge_it,edge_t>{
public:

  typedef Contour<point_it,edge_it,edge_t> base_t;
//  using base_t::points_begin;
//  using base_t::points_end;
  virtual point_it points_begin() const =0;
  virtual point_it points_end() const =0;


  bool TestPoint(const Vector_2 p) const;
};

// Contour_N specialization for 3D case
template<class point_it, class edge_it, class edge_t>
class Contour_N<point_it,Vector_3,edge_it,edge_t>: public Contour<point_it,edge_it,edge_t>{
public:

  typedef normplanes_it<edge_it> plane_it;
  typedef Contour<point_it,edge_it,edge_t> base_t;
//  using base_t::points_begin;
//  using base_t::points_end;
  virtual point_it points_begin() const =0;
  virtual point_it points_end() const =0;
  using base_t::edges_begin;
  using base_t::edges_end;

  // tests whether  a segment crosses the contour 
  // plane inside the contour and returns the crossing point (if cross!=NULL)
  bool TestEdge(const Vector_3 &v1, const Vector_3 &v2, Vector_3 *cross=NULL) const;

  /// get normal planes begin
  plane_it planes_begin() const {
    Plane_3 cpl=GetPlane();
    Vector_3 n=cpl.get_normal();
    return plane_it(n,edges_begin());
  }
  /// get normal planes end
  plane_it planes_end() const {
    return plane_it(Vector_3(0),edges_end());
  }

  /// gets the plane of the contour (takes 3 first points only)
  Plane_3 GetPlane() const {
    point_it it=points_begin();
    return Plane_3(GetNormVect(),*it);
  }

  /// gets normal vector (takes 3 first points only)
  Vector_3 GetNormVect() const;
};

// iterator at points
template<class vector_t>
class arr_point_it{
  int c; // current index
  vector_t *points; // pointer at points array
public:
  arr_point_it(int si,vector_t *spoints):c(si),points(spoints){}
  arr_point_it(const arr_point_it& other):c(other.c),points(other.points){}
   
  typedef vector_t value_type;
  
  vector_t &operator*() const {
    return points[c];
  }
  
  arr_point_it &operator++(){ // prefix
    c++;
    return *this;
  }

  arr_point_it operator++(int){
    arr_point_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const arr_point_it &other) const{
    return c!=other.c;
  }
};

template <int N=3>
class PtrContour: public Contour_N<arr_point_it<Vector_Nt<vec_type,N> > >{

public:

  typedef arr_point_it<Vector_Nt<vec_type,N> > point_it;
//  typedef Contour<point_it, generic_edge_it<point_it>, typename generic_edge_it<point_it>::edge> base_t;
  typedef Contour<point_it> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
//  using typename base_t::vector_t; 
 // using base_t::point_iterator;

  mngptr<vector_t> points; // pointer at contour vertices
  int num; // number of contour vertices

  PtrContour(mngarg<vector_t> points_=NULL, int num_=0):points(points_),num(num_){}

  point_iterator points_begin() const {
    return arr_point_it<vector_t>(0,(vector_t *)points.ptr());
  }
  
  point_iterator points_end() const {
    return arr_point_it<vector_t>(num,(vector_t *)points.ptr());
  }

  /// no checks
  vector_t GetPoint(int ind) const {
    return points[ind];
  }

  /// no checks
  int SetPoint(int i,const vector_t &p){
    points[i]=p;
    return 1;
  }

  int GetNPoints() const {
    return num;
  }
};

template <int N,int num>
class ArrayContour: public Contour_N<arr_point_it<Vector_Nt<vec_type,N> > >{

public:

  typedef arr_point_it<Vector_Nt<vec_type,N> > point_it;
//  typedef Contour<point_it, generic_edge_it<point_it>, typename generic_edge_it<point_it>::edge> base_t;
  typedef Contour<point_it> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
//  using typename base_t::vector_t; 
 // using base_t::point_iterator;
 
  vector_t points[num];

  point_iterator points_begin() const {
    return arr_point_it<vector_t>(0,(vector_t *)points);
  }
  
  point_iterator points_end() const {
    return arr_point_it<vector_t>(num,(vector_t *)points);
  }

  /// no checks
  vector_t GetPoint(int ind) const {
    return points[ind];
  }

  /// no checks
  int SetPoint(int i,const vector_t &p){
    points[i]=p;
    return 1;
  }

  int GetNPoints() const {
    return num;
  }
};

template<int N=3>
class VecContour: public Contour_N<typename vector<Vector_Nt<vec_type,N> >::const_iterator>{

public:

  typedef Contour<typename vector<Vector_Nt<vec_type,N> >::const_iterator> base_t;
  typedef typename base_t::vector_t vector_t;
  typedef typename base_t::point_iterator point_iterator;
  
  vector<vector_t> points; // vertices of the contour

  VecContour(int n=0):points(n){}

  VecContour(const vector<vector_t> &pvec):points(pvec){}  

  /// constructs a contour from input point iterator, applies additional transform p=>p*a+shift
  template<class InpIt>
  VecContour(InpIt pbeg, InpIt pend, vec_type a=1.,const vector_t &shift=vector_t()){
    for(;pbeg!=pend;++pbeg){
      points.push_back((*pbeg)*a+shift);
    }
  }

  void add(const vector_t &vect) {
    points.push_back(vect);
  }

  point_iterator points_begin() const {
    return points.begin();
  }
  point_iterator points_end() const {
    return points.end();
  }

  /// sets the ith point, if i is greated than current size
  /// appends the vector with the required amount of default Vector_3
  int SetPoint(int i,const vector_t &p){
    if(i<0)return -1;
    int n=(int)points.size();
    while(i>n){
      points.push_back(vector_t());
      n++;
    }
    if(i==n){
      points.push_back(p);
    }
    else{
      points[i]=p;
    }
    return 1;
  }

  /// no checks
  vector_t GetPoint(int ind) const {
    return points[ind];
  }

  int GetNPoints() const {
    return (int)points.size();
  }

  vector<vector_t> &GetPoints(){
    return points;
  }

  void Clear(){
    points.clear();
  }

  operator PtrContour<N>()const{
    return PtrContour<N>((vector_t *)&points[0],GetNPoints());
  }
};

/*template<class contour_t>
PtrContour<contour_t::dimension> make_PtrContour(const contour_t &cnt){
  int n=cnt.GetNPoints();
  contour_t::vector_t *points = n ? NULL : new contour_t::vector_t[n];
  n=0;
  for(contour_t::point_iterator it=cnt.points_begin(),e=cnt.points_end();it!=e;++it)
    points[n++]=*it;
  return PtrContour<contour_t::dimension>(make_mngarg(points),n);
}

template<int N,int num>
PtrContour<N> make_PtrContour(const ArrayContour<N,num> &cnt){
  return PtrContour<N>(cnt.points,num);
}*/

// calculate 2 unit vectors x, y perpeindicular to unit vector z
int set_perpendiculars(const Vector_3 &z, Vector_3 &x, Vector_3 &y);

inline Vector_2 Vector_3to2(const Vector_3 &v3, const Vector_3 &origin, const Vector_3 &x, const Vector_3 &y){
  Vector_3 tmp=v3-origin;
  return Vector_2(tmp*x, tmp*y);
}

inline Vector_3 Vector_2to3(const Vector_2 &v2, const Vector_3 &origin, const Vector_3 &x, const Vector_3 &y){
  return origin+v2[0]*x+v2[1]*y;
}

// project 3D contour to 2D using origin, x and y axis
template<class point_it>
class contour_projection_it{
  point_it it;
  Vector_3 o,x,y;
public:
  contour_projection_it(){}
  contour_projection_it(point_it it_, const Vector_3 &o_, const Vector_3 &x_, const Vector_3 &y_):it(it_),o(o_),x(x_),y(y_){}

  typedef Vector_2 value_type;
   
  Vector_2 operator*() const {
    return Vector_3to2(*it,o,x,y);
  }
  
  contour_projection_it &operator++(){ // prefix
    ++it;
    return *this;
  }

  contour_projection_it operator++(int){
    contour_projection_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const contour_projection_it &other) const{
    return it!=other.it;
  }
};

template<class contour_t>
class Contour_3to2: public Contour_N<contour_projection_it<typename contour_t::point_iterator> >{

public:

  typedef Contour<contour_projection_it<typename contour_t::point_iterator> > base_t;
  typedef typename base_t::point_iterator point_iterator;

  point_iterator b,e;

  Contour_3to2(){}
  Contour_3to2(const contour_t &cnt,const Vector_3 &o,const Vector_3 &x,const Vector_3 &y):
  b(cnt.points_begin(),o,x,y),e(cnt.points_end(),o,x,y){}

  point_iterator points_begin() const {
    return b;
  }
  
  point_iterator points_end() const {
    return e;
  }
};

///\en record to the contour cnt section of the plane plane0 by polyhedron, 
/// formed by planes iterated from iterator beg till end
///\ru записывает в контур cnt сечение плоскости plane0 многогранником, образованным плоскостями, по которым проходят итераторы beg, end
template<class plane_it, class contour_t>
int ProjectSimplex(const Plane_3& plane0, plane_it beg, plane_it end, contour_t &cnt, int nplanes=-1);

/// gets all faces of the polyhedral region defined by planes from [beg, end) of type plane_it
/// puts the faces as planar contours into container cont using push_back(...)
/// @return the number of faces
template<class out_cont_t, class plane_it>
int GetFaces(out_cont_t &cont,plane_it beg,plane_it end);

#endif
