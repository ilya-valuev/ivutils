/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.5 $
 *   $Date: 2012/06/29 10:50:13 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/region.cpp,v 1.5 2012/06/29 10:50:13 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/region.cpp,v $
$Revision: 1.5 $
$Author: valuev $
$Date: 2012/06/29 10:50:13 $
*/
/*s****************************************************************************
 * $Log: region.cpp,v $
 * Revision 1.5  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.57  2012/04/12 20:10:56  lesha
 * documentation
 *
 * Revision 1.56  2012/03/21 17:21:54  lesha
 * documentation
 *
 * Revision 1.55  2012/03/07 09:28:47  lesha
 * *** empty log message ***
 *
 * Revision 1.54  2012/02/22 03:11:59  lesha
 * *** empty log message ***
 *
 * Revision 1.53  2012/02/22 02:52:52  lesha
 * *** empty log message ***
 *
 * Revision 1.52  2012/02/17 01:06:30  lesha
 * PtrContour is added
 *
 * Revision 1.51  2011/02/26 03:15:36  lesha
 * PARDISO is included
 *
 * Revision 1.50  2010/10/22 15:43:29  biaks
 * compiled with gcc 4.4
 *
 * Revision 1.49  2010/10/19 14:11:23  valuev
 * compiled with gcc 4
 *
 * Revision 1.48  2010/10/19 11:05:58  lesha
 * *** empty log message ***
 *
 * Revision 1.47  2010/10/19 10:59:23  lesha
 * *** empty log message ***
 *
 * Revision 1.46  2010/10/05 16:56:59  lesha
 * slight modifications
 *
 * Revision 1.45  2010/09/28 08:33:14  biaks
 * add cone
 *
 * Revision 1.44  2010/07/13 17:03:18  valuev
 * layer number definition from the translation indices
 *
 * Revision 1.43  2010/05/03 21:46:57  lesha
 * GetCylinder for confined cylinder is added
 *
 * Revision 1.42  2010/03/17 17:59:45  lesha
 * uiImpulse, uiPMLFunction are removed
 *
 * Revision 1.41  2010/01/04 12:36:10  lesha
 * SpaceRegionShift is removed
 *
 * Revision 1.40  2010/01/03 13:51:52  lesha
 * TestLine is added
 *
 * Revision 1.39  2010/01/02 22:22:58  lesha
 * StretchedRegion is added, Contour::GetCenter if fixed, 2DDumping is added etc.
 *
 * Revision 1.38  2010/01/01 20:40:24  lesha
 * Region, Box_N, Sphere_N, Basis_N, MonteCarloTestContour are added
 *
 * Revision 1.37  2009/12/01 18:46:07  lesha
 * contour is modified
 *
 * Revision 1.36  2009/11/18 20:51:58  valuev
 * fixed lattice filling
 *
 * Revision 1.35  2009/10/05 07:09:10  lesha
 * make_reciprocal_vectors
 *
 * Revision 1.34  2009/10/02 17:38:40  lesha
 * simplification
 *
 * Revision 1.33  2009/09/30 14:44:58  valuev
 * Added fill_lattice function
 *
 * Revision 1.32  2009/08/25 08:48:17  lesha
 * GetCylinder is fixed
 *
 * Revision 1.31  2009/06/21 10:38:47  lesha
 * region_2.h is added
 *
 * Revision 1.30  2009/06/11 19:03:31  lesha
 * Box::TestEdge is included
 *
 * Revision 1.29  2009/05/31 11:34:06  lesha
 * Cone is added
 *
 * Revision 1.28  2009/05/23 14:41:29  lesha
 * SpaceRegionDumper is included
 *
 * Revision 1.27  2009/05/21 10:51:45  lesha
 * Box::TestRay calls Polyhedron::TestRay
 *
 * Revision 1.26  2009/05/20 12:39:29  lesha
 * Box::TestRay is modified
 *
 * Revision 1.25  2009/05/19 20:50:39  lesha
 * Box::TestRay is added
 *
 * Revision 1.24  2009/05/19 08:28:19  lesha
 * TestRay is added
 *
 * Revision 1.23  2009/05/01 11:10:33  lesha
 * Ellipsoid::TestIntersection is added
 *
 * Revision 1.22  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.21  2008/11/15 20:26:43  lesha
 * Ellipsoid is added
 *
 * Revision 1.20  2008/09/16 08:05:53  lesha
 * bug in Polyhedron_3 constructor using is fixed
 *
 * Revision 1.19  2008/07/30 10:19:41  lesha
 * Sphere::GetBoundingBox is added
 *
 * Revision 1.18  2008/03/21 08:47:05  lesha
 * bug in ConfinedBody::TestContour is fixed
 *
 * Revision 1.17  2008/02/29 09:17:42  lesha
 * ConfinedRegion is modified (pointers -> mngptr)
 *
 * Revision 1.16  2008/01/25 18:12:41  lesha
 * template instantiation is corrected
 *
 * Revision 1.15  2008/01/24 11:21:11  lesha
 * Polyhedron_3 constructor is modified
 *
 * Revision 1.14  2008/01/23 18:15:01  lesha
 * template class instantiation is corrected
 *
 * Revision 1.13  2008/01/22 10:17:06  lesha
 * Polyhedron_3 constructor is modified
 *
 * Revision 1.12  2008/01/10 06:41:57  lesha
 * GetPlane / Plate are added
 *
 * Revision 1.11  2008/01/06 03:16:08  lesha
 * common.h is excluded. Vector_3::type->vec_type
 *
 * Revision 1.10  2008/01/05 21:08:05  lesha
 * Implementation part is moved from *.h to *.cpp
 *
 * Revision 1.9  2007/09/15 08:22:13  lesha
 * TestRegionIntersection is added
 *
 * Revision 1.8  2007/09/02 15:14:35  lesha
 * Cylinder is modified
 *
 * Revision 1.7  2007/07/05 22:10:31  lesha
 * FlatRegion, RectRegion and Cylinder are added
 *
 * Revision 1.6  2007/03/05 01:10:31  valuev
 * added biolab1 pbs startup; modified interpolation for emSourceWave
 *
 * Revision 1.5  2007/03/04 04:10:45  lesha
 * ShiftSphere is added (for lattice)
 *
 * Revision 1.4  2007/02/28 21:06:46  lesha
 * SpaceRegionShift is added
 *
 * Revision 1.3  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
/// \file \brief Non-template function definitions for  region.h

#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include <algorithm>
#include "region_2.hpp"
#include "region.hpp"
#include "mnw.h"
#include "linsysn.h"


int Box::TestLine(const vector_t &p, const vector_t &dir, 
vec_type *frac, vector_t *surfp, vector_t *surfn, vec_type epsilon) const{
  Polyhedron<Box::plane_it> poly(planes_begin(),planes_end());
  return poly.TestLine(p,dir,frac,surfp,surfn,epsilon);
}

bool TestRegionsIntersection(Sphere *reg1, Sphere *reg2) {
  Vector_3 c1=reg1->get_center(), c2=reg2->get_center();
  vec_type r1=reg1->get_radius(), r2=reg2->get_radius();
  vec_type rr=r1+r2;
  return ((c2-c1).norm2()<rr*rr);
}

Polyhedron_3 *GetPolyhedronBox(const Vector_3 &p1, const Vector_3 &p2){
  Box B(p1,p2);
  Plane_3 *pl = new Plane_3[6];
  int i=0;
  for (Box::plane_it it=B.planes_begin(),e=B.planes_end();it!=e;++it){
    pl[i++]=*it;
  }
  return new Polyhedron_3(pl,pl+6);
}

Polyhedron_3 *GetPolyhedronPlane(const Vector_3 &n, const Vector_3 &pos){
  if(n==0)return NULL;
  Plane_3 *pl=new Plane_3[2];
  pl[0].init(n, pos);
  return new Polyhedron_3(pl, pl+1);
}

Polyhedron_3 *GetPolyhedronPlate(const Vector_3 &n, const Vector_3 &pos, vec_type width){
  if(n==0)return NULL;
  Plane_3 *pl=new Plane_3[3];
  pl[0].init(n, pos);
  pl[1].init(n, pos);
  Vector_3 n2;
  vec_type d;
  pl[1].get_coeff(n2, d);
  n2=-n2;
  d=-d+width;
  pl[1].set_coeff(n2, d);
  return new Polyhedron_3(pl, pl+2);
}

Polyhedron_3 *GetPolyhedronPrism(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base){
  VecContour<2> *cnt=base->get_contour();
  int sz=cnt->GetNPoints();

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  int sign=0;
  for (generic_edge_it<VecContour<2>::point_iterator> it(cnt->points_begin(), cnt->points_end()), e(cnt->points_end(), cnt->points_end(), 1); it!=e; ++it) {
    edge_t<VecContour<2>::point_iterator> ed=*it;
    Vector_2 a=ed.get_p1(), b=ed.get_p2();
    Vector_3 a3=O+a[0]*x+a[1]*y, b3=O+b[0]*x+b[1]*y;

    Vector_3 n=(top-O)%(b3-a3);
    if (sign==0){
      Vector_3 c=O-a3;
      sign=(n*c>0)? 1 : -1;
    }

    planes[sz++].init(sign*n, a3);
  }
  
  return new Polyhedron_3(planes, planes+sz);
}

Polyhedron_3 *GetPolyhedronPyramida(const Vector_3 &O, const Vector_3 &top, const Vector_3 &x, const Vector_3 &y, Polygon_2 *base){

  VecContour<2> *cnt=base->get_contour();
  int sz=cnt->GetNPoints();

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  int sign=0;
  for (generic_edge_it<VecContour<2>::point_iterator> it(cnt->points_begin(), cnt->points_end()), e(cnt->points_end(), cnt->points_end(), 1); it!=e; ++it) {
    edge_t<VecContour<2>::point_iterator> ed=*it;
    Vector_2 a=ed.get_p1(), b=ed.get_p2();
    Vector_3 a3=O+a[0]*x+a[1]*y, b3=O+b[0]*x+b[1]*y;

    Vector_3 n=(top-a3)%(b3-a3);
    if (sign==0){
      Vector_3 c=O-a3;
      sign=(n*c>0)? 1 : -1;
    }

    planes[sz++].init(sign*n, a3);
  }
  
  return new Polyhedron_3(planes, planes+sz);
}

Polyhedron_3 *GetConfinedPolyhedron(mngarg<Polyhedron_3> c1, mngarg<Polyhedron_3> c2){
  Polyhedron_3 *poly[2]={c1.first,c2.first};
  int sz=0;
  for(int i=0;i<2;i++){
    if(poly[i]){
      Polyhedron_3::plane_it it=poly[i]->planes_begin(), e=poly[i]->planes_end();
      for(; it!=e; ++it)
        sz++;
    }
  }

  Plane_3 *planes = new Plane_3[sz];
  sz=0;

  for(int i=0;i<2;i++){
    if(poly[i]){
      Polyhedron_3::plane_it it=poly[i]->planes_begin(), e=poly[i]->planes_end();
      for (; it!=e; ++it)
        planes[sz++]=*it;
    }
  }

  if(c1.second)
    delete c1.first;
  if(c2.second)
    delete c2.first;
  
  return new Polyhedron_3(planes, planes+sz);
}

Vector_3 random_direction(){
  vec_type fi=2*M_PI*(vec_type)rand()/(vec_type)RAND_MAX;
  vec_type tetta=acos(-1+2*(vec_type)rand()/(vec_type)RAND_MAX);
  vec_type sin_tetta=sin(tetta);
  return Vector_3(sin_tetta*cos(fi),sin_tetta*sin(fi),cos(tetta));
}

Plane_3* CreateSpherePlanes(const int &N, const vec_type &R, const Vector_3 &center){
  Plane_3 *SpherePlanes=new Plane_3[N+1];
  for (int i=0;i<N;i++) {
    Vector_3 n=R*random_direction();
    Vector_3 pos=center+n;
    SpherePlanes[i].init(-n,pos);
  }
  return SpherePlanes;
}

Cylinder<Circle> *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R){
  n.normalize();
  Vector_3 x(n[2]-n[1],n[0]-n[2],n[1]-n[0]);
  if(x==0)
    x=Vector_3(1,-1,0);
  x.normalize();
  Vector_3 y=n%x;
  return new Cylinder<Circle>(origin,n,x,y,make_mngarg(new Circle(R,Vector_2())));
}

ConfinedRegion<Cylinder<Circle> > *GetCylinder(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type height){
  return new ConfinedRegion<Cylinder<Circle> >(
    make_mngarg(GetCylinder(origin,n,R)),
    make_mngarg(GetPolyhedronPlate(n,origin,height))
  );
}

ConfinedRegion<Cone<Circle> > *GetCone(const Vector_3 &origin, Vector_3 n, vec_type R, vec_type h, Vector_3 npl){
  n.normalize();
  Vector_3 x(n[2]-n[1],n[0]-n[2],n[1]-n[0]);
  if(x==0)
    x=Vector_3(1,-1,0);
  x.normalize();
  Vector_3 y=n%x;

  return new ConfinedRegion<Cone<Circle> >(
    make_mngarg(new Cone<Circle>(origin,n,x,y,make_mngarg(new Circle(R,Vector_2())),h)),
    make_mngarg(GetPolyhedronPlane(npl==0 ? n : npl,origin)));
}



bool get_first_corner(const Box &box,const Plane_3 &plane,const Vector_3 &k,Vector_3 *start){
  VecContour<> cnt;
  if(ProjectSimplex(plane, box.planes_begin(), box.planes_end(), cnt)<0)
    return false;
  *start=get_min_point(k,cnt.points_begin(),cnt.points_end());
  return true;
}

bool get_first_corner(const Box &box,const Vector_3 &k,Vector_3 *start){
  if(!box.valid())
    return false;
  *start=get_min_point(k,box.points_begin(),box.points_end());
  return true;
}

bool find_cross_segment(const SpaceRegion &reg, const Vector_3& orig,const Vector_3& dir,vec_type &t1, vec_type &t2){
  Vector_3 tdir=dir;
  if(reg.TestPoint(orig)){  // the point is inside
    t2=reg.TestRay(orig,tdir);
    if(t2<0) //strange, no intersection 
      return false; 
    t1=reg.TestRay(orig,-tdir);
    if(t1<0) //strange, no intersection 
      return false;
    t1=-t1;
    return true;
  }
  int inv1=1, inv2=1;
  t1=reg.TestRay(orig,tdir);
  if(t1<0){
    t1=reg.TestRay(orig,-tdir);
    if(t1<0)
      return false;
    tdir=-tdir;
    inv1=-1;
  }
  Vector_3 p1=orig+tdir*t1;
  Vector_3 p2=p1+tdir; // test point
  if(reg.TestPoint(p2)) // test point is inside, must go in the same direction
    t2=reg.TestRay(p2,tdir);
  else{ // opposite direction
    t2=reg.TestRay(p2,-tdir);
    inv2=-1;
  }
  if(t2<0)
    t2=t1;
  else
    t2=t1+1+inv2*t2;
 
  if(inv1<0){
    t1=-t1;
    t2=-t2;
    swap(t1,t2);
  }
  return true;
}

/// Fills the space inside convex region by the lattice given by 3 elementary translations (cell) starting from origin
/// (orig). Puts the filled positions into points by push_back operation (no clearing is performed).
/// Optionally (if ipoints is not NULL) fills integer lattice indices indicating the translation numbers for each point.
/// @return the number of points filled (on success)
///         -1 no intersection of the lattice with the region found
///         -2 the region is unbounded
int fill_lattice(const SpaceRegion &reg, const Vector_3& orig,const Basis_3 &cell,vector<Vector_3> &points, vector<iVector_3> *ipoints){
  Vector_3 v1, v2, cnt;
  cnt=reg.GetBoundingBox(&v1,&v2);// always starting from region center
  if(v1.infinite() || v2.infinite())
    return -2; // can't fill unbounded region
  Vector_3 sh=cell.inv(cnt-orig);
  
  /*Box bb0(v1,v2);
  RegDumper *box_dmp=GetRegionDumper(&bb0);
  vector<VecContour> dmppoints;
  box_dmp->Dump(dmppoints);
  FILE *f1=fopen("lim.pol","wt");
  DumpVecContours(f1,dmppoints);
  fclose(f1);
  delete box_dmp;*/

  // finding limiting projection
  Box bb(v1-orig,v2-orig);
  int jmax=numeric_limits<int>::min();
  int jmin=numeric_limits<int>::max(); 
  int imin=jmin, imax=jmax;
  for(Box::point_it pit=bb.points_begin();pit!=bb.points_end();++pit){
    v1=cell.inv(*pit);
    imax=max((int)ceil(v1[0]),imax);
    imin=min((int)floor(v1[0]),imin);
    jmax=max((int)ceil(v1[1]),jmax);
    jmin=min((int)floor(v1[1]),jmin);
  } 

  /*jmax=20;
  jmin=-20;
  imax=20;
  imin=-20;*/
  /*f1=fopen("tst.d","wt");
  FILE *f2=fopen("tstc.d","wt");
  FILE *f3=fopen("tstn.d","wt");*/


  int i[2]={(int)floor(sh[0]),(int)ceil(sh[0])};
  if(i[0]==i[1])
    i[1]++;
  for(int id=0;id<2;id++){
    int ncrossi=0;
    do{
      int ncrossj=0;
      int j[2]={(int)floor(sh[1]),(int)ceil(sh[1])};
      if(j[0]==j[1])
        j[1]++;
      for(int jd=0;jd<2;jd++){
        do{
          Vector_3 cntj=orig+i[id]*cell[0]+j[jd]*cell[1]+sh[2]*cell[2];
          vec_type t1,t2;
          if(!find_cross_segment(reg,cntj,cell[2],t1,t2)){ //no intersection
            if(ncrossj)
              break;
          }
          else{
            ncrossj++;
            
            /*Vector_3 v1=orig+i[id]*cell[0]+j[jd]*cell[1]+(t1+sh[2])*cell[2];
            fprintf(f2,"%g %g %g \"%d %d\"\n",v1[0],v1[1], v1[2],i[id],j[jd]);
            v1=orig+i[id]*cell[0]+j[jd]*cell[1]+(t2+sh[2])*cell[2];
            fprintf(f2,"%g %g %g \"%d %d\"\n\n\n",v1[0],v1[1], v1[2],i[id],j[jd]);
            fflush(f2);*/

            int ks=(int)ceil(t1+sh[2]), ke=(int)floor(t2+sh[2]);
            for(int k=ks;k<=ke;k++){
              Vector_3 atompos=orig+i[id]*cell[0]+j[jd]*cell[1]+k*cell[2];
              points.push_back(atompos);
              if(ipoints)
                ipoints->push_back(iVector_3(i[id],j[jd],k));
            }
          }
          if(jd){
            j[jd]++;
            if(j[jd]>jmax)
              break;
          }
          else{
            j[jd]--;
            if(j[jd]<jmin)
              break;
          }
        }while(1);
      }// jd
      if(ncrossj==0 && ncrossi!=0)
        break;  // went out of intersection region
      ncrossi+=ncrossj;
      if(id){
        i[id]++;
        if(i[id]>imax)
          break;
      }
      else{
        i[id]--;
        if(i[id]<imin)
          break;
      }

    }while(1);
  }
  //fclose(f1);
  return (int)points.size();
}

class CompareNorms{
public:
  bool operator()(const Vector_3 a, const Vector_3 b){
    return a.norm2()<b.norm2();
  }
};

Vector_3 *make_reciprocal_vectors(const vec_type a1, const vec_type a2, const vec_type max, int &num){

  vec_type kmax=1/max;

  vec_type a[2]={a1,a2};
  int n[2];
  for(int i=0;i<2;i++){
    a[i]=1/a[i];
    n[i]=int(floor(kmax/a[i]));
  }
  num=(2*n[0]+1)*(2*n[1]+1);
  Vector_3 *v=new Vector_3[num];
  int k=0;
  for(int i=-n[0];i<=n[0];i++){
    for(int j=-n[1];j<=n[1];j++){
      v[k++]=Vector_3(a[0],0,0)*i+Vector_3(0,a[1],0)*j;
    }
  }
  sort(v,v+num,CompareNorms());

  return v;
}

template class Box_N<2>;
template class Box_N<3>;
template class Sphere_N<2>;
template class Sphere_N<3>;

template class Polyhedron<Box::plane_it>;
template class Polyhedron<Plane_3 *>;

template class Cylinder<Circle>;
template class Cone<Circle>;
template class StretchedRegion<SpaceRegion>;
template class ConfinedRegion<SpaceRegion>;

template vec_type Box::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type Sphere::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type Polyhedron<Box::plane_it>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type Polyhedron<Plane_3 *>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;


template vec_type Circle::TestContour(const PtrContour<2> &cnt, VecContour<2> *subcont, Vector_2 *subcenter) const;

template vec_type Cylinder<Circle>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
//template vec_type Cone<Circle>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;
template vec_type StretchedRegion<Region<3> >::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Sphere>;
template vec_type ConfinedRegion<Sphere>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Cylinder<Circle> >;
template vec_type ConfinedRegion<Cylinder<Circle> >::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template class ConfinedRegion<Cone<Circle> >;

template class ConfinedRegion<MNW>;
template vec_type ConfinedRegion<MNW>::TestContour(const PtrContour<3> &, VecContour<> *, Vector_3 *) const;

template int MonteCarloTestLine(const Region<3> *reg, const Vector_Nt<vec_type,3> &p, const Vector_Nt<vec_type,3> &dir, vec_type *frac, 
Vector_Nt<vec_type,3> *surfp, Vector_Nt<vec_type,3> *surfn, vec_type length, int nt);
template vec_type MonteCarloTestContour(const Region<2> *reg, const VecContour<2> &cnt, Vector_Nt<vec_type,2> *subcenter, int nt);
template vec_type MonteCarloTestContour(const Region<3> *reg, const VecContour<3> &cnt, Vector_Nt<vec_type,3> *subcenter, int nt);

template int MakeBoxDomains(const Box &box,int nproc, int nx, int ny, int nz, vector<Box>::iterator result);
