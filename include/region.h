/*s****************************************************************************
 * $Log: region.h,v $
 * Revision 1.6  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.99  2012/06/07 17:22:23  valuev
 * started adding scalable dipole
 *
 * Revision 1.98  2012/05/13 15:12:22  lesha
 * russian documentation reverted back
 *
 * Revision 1.97  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.96  2012/03/21 20:59:08  lesha
 * shift polyhedron is fixed
 *
 * Revision 1.95  2012/03/21 17:16:54  lesha
 * documentation
 *
 * Revision 1.94  2012/02/17 00:17:36  lesha
 * PtrContour is added
 *
 * Revision 1.93  2012/01/10 01:07:11  lesha
 * SurfProject is added
 *
 * Revision 1.92  2011/11/29 23:27:51  lesha
 * Dump of MNW is added
 *
 * Revision 1.91  2011/09/02 03:16:25  lesha
 * Polyhedron_3 constructor bug is fixed
 *
 * Revision 1.90  2011/01/26 00:59:12  lesha
 * TestLine is added to Inverse and StretchedRegion
 *
 * Revision 1.89  2010/10/21 16:06:46  biaks
 * add brief tag and comments
 *
 * Revision 1.88  2010/10/15 20:55:03  lesha
 * *** empty log message ***
 *
 * Revision 1.87  2010/10/15 20:08:52  lesha
 * *** empty log message ***
 *
 * Revision 1.86  2010/10/15 16:19:38  lesha
 * *** empty log message ***
 *
 * Revision 1.85  2010/10/15 14:36:16  lesha
 * *** empty log message ***
 *
 * Revision 1.84  2010/10/06 11:28:16  biaks
 * add comments
 *
*******************************************************************************/
#ifndef _REGION_H
#define _REGION_H

/// \en @file region.h \brief A collection of objects describing various 3D geometrical bodies.
/// \ru @file region.h \brief Набор объектов, описывающих различные 3D геометрические тела.


#include "region_2.h"


// polyhedron traits
template<class poly_t>
struct poly_traits{
  typedef void plane_it; 
};

/// orthogonal box, see documentation to Box_N
class Box: public Box_N<3>{
  
public:
  typedef polyhedron_cat category;

  Box(){}

  Box(const Vector_3 &sp1, const Vector_3 &sp2){
    init(sp1,sp2);
  }

  Box(vec_type x1, vec_type x2,vec_type y1, vec_type y2,vec_type z1, vec_type z2){
    init(Vector_3(x1,y1,z1),Vector_3(x2,y2,z2));
  }

  ///\en Returns a set of bit flags specifying which functions
  ///    are implemented for this region (see \ref RegFlags).
  ///    \return HAS_TESTPOINT|HAS_TESTLINE|HAS_TESTCONTOUR|HAS_CLONE|HAS_VOLUME
  virtual int GetFlags() const {
    return HAS_TESTPOINT|HAS_TESTLINE|HAS_TESTCONTOUR|HAS_CLONE|HAS_VOLUME;
  }

  class plane_it {
    friend class Box;
    const Box *parent;
    int i;
    plane_it(const Box *sparent, int si=0):parent(sparent), i(si){}
  public:
    typedef Plane_3 value_type;
    /// default constructor pointing to the end
    plane_it():parent(NULL),i(6){}
    
    plane_it(const plane_it &other):parent(other.parent),i(other.i){}
    plane_it &operator++(){
      i++;
      return *this;
    }
    plane_it operator++(int){
      plane_it tmp=*this;
      ++*this;
      return tmp;
    }
    /// plane construction
    Plane_3 operator*() const{
      Vector_3 n;
      n[i%3]= (i>2? -1: 1); // other coords =0
      return Plane_3(n,(i>2? parent->p2 : parent->p1));
    }
    bool operator!=(const plane_it &other) const{
      return i!=other.i;
    }
  };

  plane_it planes_begin() const{
    return plane_it(this,0);
  }
  plane_it planes_end() const{
    return plane_it(this,6);
  }

  virtual vec_type Volume() const{
    vec_type vol=1;
    for(int i=0;i<3;i++)vol*=sz[i];
    return vol;
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  ///\en Returns the area fraction of the contour part that is inside the region 
  /// and subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

// polyhedron traits
template<>
struct poly_traits<Box>{
  typedef Box::plane_it plane_it; 
};

#if 0
/// parallelepiped
class ShearBox: public ShearBox_N<3>{
public:
  typedef polyhedron_cat category;

  ShearBox(){}

  /// constructs a shear box expressed via a Box in the own basis
  ShearBox(const Vector_3 &sp1, const Vector_3 &sp2, const Basis_3 &basis=Basis_3()){
    this->basis=basis;
    init(sp1,sp2);
  }

  /// constructs from origin and 3 sides given in the orthogonal basis
  ShearBox(const Vector_3 &origin, const Vector_3 &v0, const Vector_3 &v1, const Vector_3 &v2){
    Vector_3 p=v0+v1+v2;
    basis=Basis_3(v0.unit(),v1.unit(),v2.unit());
    init(basis.inv(origin),basis.inv(p));
  }

  class plane_it {
    friend class ShearBox;
    const ShearBox *parent;
    int i;
    /// end=6
    plane_it(const ShearBox *sparent, int si=0):parent(sparent), i(si){}
  public:
    typedef Plane_3 value_type;
    /// default constructor pointing to the end
    plane_it():parent(NULL),i(6){}
    
    plane_it(const plane_it &other):parent(other.parent),i(other.i){}
    plane_it &operator++(){
      i++;
      return *this;
    }
    plane_it operator++(int){
      plane_it tmp=*this;
      ++*this;
      return tmp;
    }
    /// plane construction
    Plane_3 operator*() const{
      Vector_3 n;
      n[i%3]= (i>2? -1: 1); // other coords =0
      Vector_3 p=(i>2? parent->p2 : parent->p1);
      return Plane_3(parent->basis.cov(n),n*p);
    }
    bool operator!=(const plane_it &other) const{
      return i!=other.i;
    }
  };

  plane_it planes_begin() const{
    return plane_it(this,0);
  }
  plane_it planes_end() const{
    return plane_it(this,6);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};
#endif

/// plyhedron defined as a set of planes
template <class plane_itt>
class Polyhedron: public SpaceRegion{
protected:
  plane_itt b, e; // start and end iterator at plane sequence which confine the polyhedron
  Box box;
public:
  typedef polyhedron_cat category;
  typedef plane_itt plane_it;

  Polyhedron() {};

  Polyhedron(plane_it beg, plane_it end):box(Vector_3(-VEC_INFTY),Vector_3(VEC_INFTY)){
    init(beg,end);
  }

  plane_it planes_begin() const{
    return b;
  }

  plane_it planes_end() const{
    return e;
  }

  // start and end iterator at plane sequence which confine the polyhedron
  void init(plane_it beg, plane_it end);

  /// returns true if the point is in the the region
  /// plane normals are pointing inside the region!!!!
  virtual bool TestPoint(const Vector_3 &p) const;

  /// gets the distance between a point and the closest plane
  /// (which is the distance to the polyhedron surface for internal points)
  vec_type MinPlaneDist(const Vector_3 &pos, plane_it *mit=NULL);
  
  /// returns the area fraction of the contour part that is inside the region 
  /// and subcontour which is inside the region (if pointer is not NULL)
  /// and the center of the subcontour (if pointer is not NULL)
  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  // find polyhedron surface point surfp which is closest to the given point p and normal to polyhedron surfn at surfp
  // return distance between p and surfp
  virtual vec_type SurfProject(const vector_t &p, vector_t *surfp=NULL, vector_t *surfn=NULL) const;

  virtual Vector_3 GetBoundingBox(Vector_3 *v1, Vector_3 *v2) const {
    return box.GetBoundingBox(v1, v2);
  }

  /// cloning is not possible for arbitrary plane iterator
  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const {
    return NULL;
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

// polyhedron traits
template <class plane_itt>
struct poly_traits<Polyhedron<plane_itt> >{
  typedef typename Polyhedron<plane_itt>::plane_it plane_it; 
};


template<>
struct value_type_cast<Plane_3 *>{
  typedef Plane_3 value_type;
};

/// polyhedron as array of planes
class Polyhedron_3: public Polyhedron<Plane_3 *>{
protected:
  mngptr<Plane_3> pptr;
public:
  typedef polyhedron_cat category;
  typedef Polyhedron<Plane_3 *> base_t;
  typedef base_t::plane_it plane_it;

  Polyhedron_3(){}

  /// copy constructor
  Polyhedron_3(const Polyhedron_3 &other){
    size_t sz=other.e-other.b;
    if(!sz)
      return;
    pptr.reset(new Plane_3[sz],0x8);
    for(size_t i=0;i<sz;i++)
      pptr.ptr()[i]=other.b[i];

    Polyhedron<Plane_3 *>::init(pptr.ptr(),pptr.ptr()+sz);
    box=other.box;
  }

  Polyhedron_3(Plane_3 *beg, Plane_3 *end, bool managed=true):base_t(beg,end),pptr(beg,managed ? 0x8 : 0){}

  /// the pointers must be valid for the whole lifetime of the class
  /// beg and end start and end pointers at array of planes
  /// managed: false if this array is not managed, true if this array will be delete if not used anymore
  void init(Plane_3 *beg, Plane_3 *end, bool managed=true){
    base_t::init(beg,end);
    pptr.reset(beg,managed ? 0x8 : 0);
  }
  
  /// this constructor copies the planes to the current class
  /// the managed flag is ignored (as if were always true)
  template<class plane_itt>
  void init(plane_itt beg, plane_itt end, bool managed){
    plane_itt it;
    size_t sz=0;
    for(it=beg;it!=end;++it)sz++;
    Plane_3 *data= sz? new Plane_3[sz] : NULL; // incomplete deleting
      
    pptr.reset(data, sz? 0x8: 0);
    sz=0;
    // copying
    for(it=beg;it!=end;++it){
      data[sz]=*it;
      sz++;
    }
    Polyhedron<Plane_3 *>::init(data,data+sz);
  }

  /// this constructor copies the planes to the current class
  /// the managed flag is ignored (as if were always true)
  template<class plane_itt>
  Polyhedron_3(plane_itt beg, plane_itt end, bool managed){
    init(beg,end,managed);
  }

  Polyhedron_3(const vector<Plane_3> &planes){
    init(planes.begin(),planes.end(),true);
  }

  Polyhedron_3(const Box &b){
    init(b.planes_begin(),b.planes_end(),true);
  }

  Polyhedron_3(const Plane_3 &pl){
    init(&pl,&pl+1,true);
  }

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  /// @return a copy of itself optionally shifted in space with the shift vector
  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const {
    size_t sz=e-b;
    if(!sz)
      return NULL;
    Plane_3 *nplanes= new Plane_3[sz];
    for(size_t i=0;i<sz;i++){
      nplanes[i]=b[i];
      nplanes[i].shift(shift);
    }
    return new Polyhedron_3(nplanes,nplanes+sz);
  }
};

// polyhedron traits
template<>
struct poly_traits<Polyhedron_3>{
  typedef Polyhedron_3::plane_it plane_it; 
};

// 3D sphere
class Sphere: public Sphere_N<3>{

public:
  typedef sphere_cat category;

  Sphere(vec_type R, const Vector_3 &center): Sphere_N<3>(R,center){}

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual Region<3> *Clone(const vector_t& shift=vector_t()) const{
    return new Sphere(R,center+shift);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en
/// Cylinder with arbitrary 2D base of the type base_t
/// Cylinder axis passes the point origin and directed along vector n
/// Axis x and of 2D plane which contains the base are directed along vector x and y
/// Vector n is perpendicular to x and y (axis is perpendicular to the base)
/// "Tilted" cylinder with axis not perpendicular to the base,
/// can be specified using template StretchedRegion
///\ru
/// Цилиндр с произвольным двумерным основанием типа base_t
/// Ось цилиндра прохожит через точку origin и направлена вдоль оси n
/// Оси x и y двумерной плоскости, на которой лежит основание, направлены вдоль векторов x и y
/// Подразумевается, что вектор n перпендикулярен x и y (основание перпендикулярно оси циллиндра)
/// "Косой" циллиндр, у которого заданное основание не перпендикулярно оси циллиндра, 
/// нужно задавать с помощью шаблона StretchedRegion
template<class base_tt>
class Cylinder: public SpaceRegion{
public:
  typedef base_tt base_t;

  Vector_3 origin, n, x, y;
  ///\ru  основание
  mngptr<base_t> base; ///<\en base

  ///\en project point v at base 2D coordinate system and record to height distance to this base plane
  ///\ru проецирует точку на плоскость основания и записывает в height высоту над этой плоскостью
  Vector_2 Vector_3to2(const Vector_3 &v3, vec_type *height=NULL) const{
    Vector_3 tmp=v3-origin;
    if(height)*height=tmp*n;
    return Vector_2(tmp*x, tmp*y);
  }

  
  ///\en move point from 2D base coordinate system to 3D space and shift from the base plane at distance height
  ///\ru переводит точку обратно в трехмерное пространство и поднимает на высоту height
  Vector_3 Vector_2to3(const Vector_2 &v2, const vec_type height=0) const{
    return origin+v2[0]*x+v2[1]*y+height*n;
  }

  Cylinder(const Vector_3 &origin, const Vector_3 &n, const Vector_3 &x, const Vector_3 &y, mngarg<base_t> base){
    init(origin, n, x, y, base);
  }

//  Cylinder(Cylinder &other): origin(other.origin), n(other.n), x(other.x), y(other.y) {
//    base.reset(make_mngarg(new base_t(*(other.base))));
//  }

  // n, x, y must be orthonormal
  void init(const Vector_3 &origin_, const Vector_3 &n_, const Vector_3 &x_, const Vector_3 &y_, mngarg<base_t> base_){
    origin=origin_, n=n_, x=x_, y=y_, base.reset(base_);
  }

  base_t *get_base() const{
    return base.ptr();
  }

  bool TestPoint(const Vector_3 &p) const{
    return base->TestPoint(Vector_3to2(p));
  }

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const {
    return new Cylinder(origin+shift,n,x,y, make_mngarg(new base_t(*base)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

// Cone with arbitrary 2D base of the type base_tt
// implementation of functions is not tested
template<class base_tt>
class Cone: public Cylinder<base_tt>{
public:  
  typedef base_tt base_t;
  typedef typename Region<Cylinder<base_tt>::dimension>::vector_t vector_t;
public:
  vec_type L; //distance between origin and top in n direction 
  using Cylinder<base_tt>::origin;
  using Cylinder<base_tt>::n;
  using Cylinder<base_tt>::x;
  using Cylinder<base_tt>::y;
  using Cylinder<base_tt>::base;
public:

  Cone(const Vector_3 &origin, const Vector_3 &n, const Vector_3 &x, const Vector_3 &y, mngarg<base_t> base, vec_type L_):
    Cylinder<base_t>(origin,n,x,y,base),L(L_){}

  bool TestPoint(const Vector_3 &p) const{
    vec_type h;
    Vector_2 p2=Cylinder<base_t>::Vector_3to2(p,&h);
    if(h>=L)return false; // behind top
    Vector_2 o2=Cylinder<base_t>::Vector_3to2(origin);
    Vector_2 dist=p2-o2;
    dist*=L/(L-h);
    return base->TestPoint(o2+dist);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;

  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const {
    return new Cone(origin+shift,n,x,y, make_mngarg(new base_t(*base)),L);
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en Region which is intersection of polyhedrin with arbitrary other region.
/// Can be used to specify hemisphere, finite cylinders etc.
///\ru тело, получающееся в результате пересечения многогранника с каким-нибудь другим произвольным телом
/// может использоваться для получения полусфер, ограниченных циллиндров и т.д.
template<class reg_tt>
class ConfinedRegion: public SpaceRegion{
public:
  typedef reg_tt reg_t;

  mngptr<reg_t> reg; // arbitrary region
  mngptr<Polyhedron_3> poly; // confined polyhedron

  ConfinedRegion(mngarg<reg_t> reg_, mngarg<Polyhedron_3> poly_): reg(reg_), poly(poly_){}

  bool TestPoint(const Vector_3 &p) const{
    return reg->TestPoint(p) && poly->TestPoint(p);
  }

  template <class contour_t>
  vec_type TestContour(const contour_t &cnt, VecContour<> *subcont=NULL, Vector_3 *subcenter=NULL) const;

  virtual vec_type TestPtrContour(const PtrContour<3> &cnt, VecContour<3> *subcont=NULL, vector_t *subcenter=NULL) const{
    return TestContour(cnt,subcont,subcenter);
  }

  virtual int TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
  vector_t *surfp=NULL, vector_t *surfn=NULL, vec_type epsilon=0) const;


  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const {
    return new ConfinedRegion(make_mngarg((reg_t *)reg->Clone(shift)),make_mngarg((Polyhedron_3 *)poly->Clone(shift)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

///\en Complement to the region to the whole 3D space
///\ru Дополнение к телу в трехмерном пространстве
template<class reg_tt>
class Inverse: public SpaceRegion{
public:
  typedef reg_tt reg_t;
  mngptr<reg_t> reg;

  Inverse(mngarg<reg_t> reg_=NULL): reg(reg_){}

  bool TestPoint(const Vector_3 &p) const{
    return !(reg->TestPoint(p));
  }

  virtual int TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
  Vector_3 *surfp=NULL, Vector_3 *surfn=NULL, vec_type epsilon=0) const{
    int res=reg->TestLine(p,dir,frac,surfp,surfn,epsilon);
    if(surfn){
      for(int i=0;i<res;i++)
        surfn[i]*=-1;
    }
    return res;
  }

  virtual SpaceRegion *Clone(const Vector_3& shift=Vector_3()) const{
    if(!reg.ptr())return NULL;
    return new Inverse(make_mngarg((reg_t *)reg->Clone(shift)));
  }

  virtual RegDumper<3>* CreateDumper() const {
    return GetRegionDumper(this);
  }
};

inline bool TestRegionsIntersection(SpaceRegion *reg1, SpaceRegion *reg2){return true;}

inline bool TestRegionsIntersection(Box *reg1, Box *reg2) {
  return reg1->TestBox(*reg2);
}

bool TestRegionsIntersection(Sphere *reg1, Sphere *reg2);

///\en returns "polyhedron" formed by one plane (half-space)
///\ru возвращает "многогранник", образованный одной плоскостью (полупространство)
Polyhedron_3 *GetPolyhedronPlane(const Vector_3 &n, const Vector_3 &pos);
///\en returns "polyhedron", formed by 2 parallel planes (infinite plane)
///\ru возвращает "многогранник", образованный двумя параллельным плоскостями (бесконечная пластинка)
Polyhedron_3 *GetPolyhedronPlate(const Vector_3 &n, const Vector_3 &pos, vec_type width);
///\en returns polyhedron - regular parallelipiped
///\ru возвращает многогранник - прямоугольный параллелепипед
Polyhedron_3 *GetPolyhedronBox(const Vector_3 &p1, const Vector_3 &p2);
///\en returns polyhedron - prism
///\ru возвращает многогранник - призму
Polyhedron_3 *GetPolyhedronPrism(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base);
///\en returns polyhedron - pyramid
/// O is center of the polygon in the base, top is top of the pyramid, 
/// x and y are axis directions, polygon will be created relatively to these axis,
/// base is polygon , confinment is polyhedron which confines pyramid (if it is absent then pyramid is infinite in one direction)
///\ru возвращает многогранник - пирамиду
/// O - центр многоугольника, top - острие пирамиды, 
/// x и y - направления осей, относительно которых строится многоугольник,
/// base - многоугольник, confinment - многогранник, ограничивающий пирамиду (если его нет - пирамида бесконечна)
Polyhedron_3 *GetPolyhedronPyramida(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base);
///\en returns polyhedron which is intersection of two other polyhedra
///\ru возвращает многогранник, являющийся пересечением двух аргументов
Polyhedron_3 *GetConfinedPolyhedron(mngarg<Polyhedron_3> c1, mngarg<Polyhedron_3> c2);
///\en returns region which is intersection of some other region and polyhedron
template<class reg_t>
ConfinedRegion<reg_t> *GetConfinedRegion(mngarg<reg_t> c1,mngarg<Polyhedron_3> c2){
  return new ConfinedRegion<reg_t>(c1,c2);
}
///\en returns polyhedron for N sides which approximating sphere
///\ru возвращает многогранник из N граней, интерполирующий сферу
Plane_3* CreateSpherePlanes(const int &N, const vec_type &R, const Vector_3 &center);
///\en returns cylinder, origin is some point at cylinder axis, n is axis direction, R is radius
///\ru origin - произвольная точка на оси, n - направление оси, R - радиус циллиндра
Cylinder<Circle> *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R);
///\en height is height of finite cylinder, height is measured from origin in direction n
///\ru height - высота ограниченного цилиндра, отмеряемая от точки origin в направлении n
ConfinedRegion<Cylinder<Circle> > *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type height);
///\en height is height of cone, height is measured from origin in direction n
ConfinedRegion<Cone<Circle> > *GetCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type h, Vector_3 npl=0);

/// Domain decomposition of the Box, nproc<0 is the size request for np
template<class inp_it>
int MakeBoxDomains(const Box &box,int nproc, int nx, int ny, int nz, inp_it result);


/// gets the point having minimal projection on the direction k
/// infinite vectors are ignored
/// returns an infinite vector for void sequence
template<class point_it> 
Vector_3 get_min_point(const Vector_3 &k, point_it beg, point_it end){
  vec_type dmin=VEC_INFTY;
  Vector_3 res(VEC_INFTY);
  for(;beg!=end;++beg){
    Vector_3 v=*beg;
    if(!v.infinite()){
      vec_type d=k*v;
      if(d<dmin){
        dmin=d;
        res=v;
      }
    }
  }
  return res;
}

/// calculates the intersection contour of the plane with the box, returns true if the intersection is nonvoid
/// *start is updated with the contour point having minimal projection on k direction 
bool get_first_corner(const Box &box,const Plane_3 &plane,const Vector_3 &k,Vector_3 *start);


/// *start is updated with the box point having minimal projection on k direction 
/// @return true if box is nonvoid, false otherwise
bool get_first_corner(const Box &box,const Vector_3 &k,Vector_3 *start);

// get arbitrary unit vector in 3D
Vector_3 random_direction();

/// Find the intersection of a convex region with the line l(t)=origin+dir*t.
/// @return true and the segment coordinates in the line system into t1 and t2 (t1<=t2) or
///         false if no intersection found
bool find_cross_segment(const SpaceRegion &reg, const Vector_3& orig,const Vector_3& dir,vec_type &t1, vec_type &t2);


/// Fills the space inside convex region by the lattice given by 3 elementary translations (cell) starting from origin
/// (orig). Puts the filled positions into points by push_back operation (no clearing is performed).
/// Optionally (if ipoints is not NULL) fills integer lattice indices indicating the translation numbers for each point.
/// @return the number of points filled (on success)
///         -1 no intersection of the lattice with the region found
///         -2 the region is unbounded
int fill_lattice(const SpaceRegion &reg, const Vector_3& orig,const Basis_3 &cell,vector<Vector_3> &points, vector<iVector_3> *ipoints=NULL);

// until Basis_Nt and calculating inverse matrix are not implemented, 
// works only for 2D rectangular case
Vector_3 * make_reciprocal_vectors(const vec_type a1, const vec_type a2, const vec_type max, int &num);

#endif
