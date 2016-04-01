/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.3 $
 *   $Date: 2012/06/29 10:50:13 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/region.hpp,v 1.3 2012/06/29 10:50:13 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/region.hpp,v $
$Revision: 1.3 $
$Author: valuev $
$Date: 2012/06/29 10:50:13 $
*/
/*s****************************************************************************
 * $Log: region.hpp,v $
 * Revision 1.3  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.44  2012/04/12 20:10:56  lesha
 * documentation
 *
 * Revision 1.43  2012/03/21 17:23:21  lesha
 * documentation
 *
 * Revision 1.42  2012/01/10 01:07:12  lesha
 * SurfProject is added
 *
 * Revision 1.41  2011/09/29 07:12:02  lesha
 * ConfinedRegion::TestLine - other case is considered
 *
 * Revision 1.40  2011/03/18 03:34:14  lesha
 * TestLine is realized for ConfinedRegion (not fully)
 *
*******************************************************************************/

/// \file \brief Template function definitions for region.h

#include "region.h"
#include "utiltl.h"
#include "linsysn.h"

#if !defined(DUMP)
template<class reg_t>
RegDumper *GetRegionDumper(const reg_t *reg){
  return NULL;
}
#endif


template <class contour_t>
vec_type Box::TestContour(const contour_t &cnt, VecContour<3> *subcont, vector_t *subcenter) const{
  Polyhedron<plane_it> poly(planes_begin(),planes_end());
  return poly.TestContour(cnt,subcont,subcenter);
}

template <class plane_itt>
void Polyhedron<plane_itt>::init(plane_it beg, plane_it end){
  b=beg;
  e=end;
  /// checking for coordinate planes
  Vector_3 n, p1(-VEC_INFTY), p2(VEC_INFTY);
  vec_type d;
  plane_it it;
  int i;
  for(it=b;it!=e;it++){
    const Plane_3 &p=*it;
    p.get_coeff(n,d);
    for(i=0;i<3;i++){
      if(fabs(fabs(n[i])-1.)<VEC_ZERO){ // coordinate direction
        break;
      }
    }
    if(i<3){ // adding as bound
      if(n[i]>0){ // left bound
        if(p1[i]<-d)p1[i]=-d;
      }
      else{  // right bound
        if(p2[i]>d)p2[i]=d;
      }
    } 
  }
  box.init(p1,p2);
}

template <class plane_itt>
bool Polyhedron<plane_itt>::TestPoint(const Vector_3 &p) const{
  // checking bounding box
  if(!box.TestPoint(p))return false;
  plane_it it;
  int res=1;
  for(it=b;it!=e;it++){
    if((*it).distance(p)<0){
      res=0;
      break;
    }
  }
  if(!res)return false;
  return true;
}

template <class plane_itt>
vec_type Polyhedron<plane_itt>::MinPlaneDist(const Vector_3 &pos, plane_it *mit){
  plane_it it=planes_begin(), e=planes_end();
  vec_type md=-1.;
  for(it=b;it!=e;it++){
    vec_type d=fabs((*it).distance(pos));
    if(md<0 || d<md){
      md=d;
      if(mit)*mit=it;
    }
  }
  return md;
}

template <class plane_itt>
template <class contour_t>
vec_type Polyhedron<plane_itt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  /// testing bounding box
  Vector_3 bc1, bc2;
  cnt.GetBoundingBox(&bc1,&bc2);
  Vector_3 p1=box.get_p1(), p2=box.get_p2();
  for(int i=0;i<3;i++){
    if(bc1[i]>p2[i] || bc2[i]<p1[i])return 0.;
  }
  aggregate_it<plane_it, typename contour_t::plane_it, Plane_3> a_b(b,e,cnt.planes_begin()),a_e(e,e,cnt.planes_end());
  // solving for intersection
  VecContour<> s_cont;
  int res=ProjectSimplex(cnt.GetPlane(),a_b,a_e,s_cont);
  if(res==1){
    vec_type a=s_cont.Area();
    if(subcenter)*subcenter=s_cont.GetCenter();
    if(subcont)subcont->GetPoints().swap(s_cont.GetPoints());
    return a; 
  }
  return 0.;
}

template <class plane_itt>
int Polyhedron<plane_itt>::TestLine(const vector_t &p, const vector_t &dir, 
vec_type *frac, vector_t *surfp, vector_t *surfn, vec_type epsilon) const{

  vec_type tleft=-numeric_limits<vec_type>::max(), tright=numeric_limits<vec_type>::max();
  Vector_3 nleft, nright;
  
  for(plane_it pi=b;pi!=e && tleft<=tright ;++pi){
    Vector_3 vnorm=(*pi).get_normal();
    vec_type dprod=-dir*vnorm; // '-' takes the opposite normal direction in ivutils into account
    // dprod>0 - я вхожу в плоскость изнутри, dprod<0 - снаружи
    vec_type dist=(*pi).distance(p); // distance to the plane

    if(fabs(dprod)<VEC_ZERO){ // dir is parallel to the plane
      if(dist>=0)
        continue; // dir is all inside, OK
      else
        return 0; // dir is all outside, no solution
    }
    vec_type t=dist/dprod;
    if(dprod>0){ // right bound, searching minimum
      if(tright>t){
        tright=t;
        nright=vnorm;
      }
    }
    else{// left bound, searching maximum
      if(tleft<t){
        tleft=t;
        nleft=vnorm;
      }
    }
  }
  if(tleft>tright)
    return -1; // no solution on this line
  frac[0]=tleft;
  frac[1]=tright;

  if(surfp){
     surfp[0]=p+dir*tleft;
     surfp[1]=p+dir*tright;
  }
  if(surfn){
    surfn[0]=-nleft;
    surfn[1]=-nright;
  }
  return 1+(tleft!=tright);
}

template <class plane_itt>
vec_type Polyhedron<plane_itt>::SurfProject(const vector_t &p, vector_t *surfp, vector_t *surfn) const{

  bool in=TestPoint(p);

  vec_type t = in ? numeric_limits<vec_type>::max() : -numeric_limits<vec_type>::max();
  Vector_3 n;
  
  for(plane_it pi=b;pi!=e;++pi){
    Vector_3 vnorm=(*pi).get_normal();
    // dist>0 - я вхожу в плоскость изнутри, dist<0 - снаружи
    vec_type dist=(*pi).distance(p); // distance to the plane

    if(in && dist>=0){ // right bound, searching minimum
      if(t>dist){
        t=dist;
        n=vnorm;
      }
    }
    else if(!in && dist<0){// left bound, searching maximum
      if(t<dist){
        t=dist;
        n=vnorm;
      }
    }
  }

  if(surfp)
     *surfp=p+n*t;

  if(surfn)
    *surfn=-n;

  return t;
}

template <class contour_t>
vec_type Sphere::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{

  if(subcenter)
    *subcenter=0;

  Vector_3 nn=cnt.GetNormVect(); // normal to contour
  nn.normalize();
    
  vec_type dist=(*(cnt.points_begin())-center)*nn; // "distance" from center of sphere to contour plane (can be negative)
  vec_type r2=(R*R-dist*dist); // radius^2 of this circle forming by intersection between contour plane and sphere
  if (r2<=0) return 0; // contour plane does not intersect sphere

  Vector_3 c=center+dist*nn; // center of circle 
  Circle C(sqrt(r2),Vector_2());
  Vector_3 ax,ay;
  set_perpendiculars(nn,ax,ay);
  Contour_3to2<contour_t> cnt2(cnt,c,ax,ay);
  Vector_2 subcenter2;
  vec_type area=C.TestContour(cnt2,NULL,subcenter ? &subcenter2 : NULL);
  if(subcenter && area)
    *subcenter=Vector_2to3(subcenter2,c,ax,ay);

  return area;
}

template<class base_tt>
template<class contour_t>
vec_type Cylinder<base_tt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  if(cnt.GetNPoints()==0)
    return 0;

  Vector_3 ncnt=cnt.GetNormVect();
  ncnt.normalize();
  vec_type coef=ncnt*n;

  if(acccomp(coef,0.)){ // контур параллелен оси циллиндра, его нельзя спроектировать на основание
    Vector_3 k=ncnt%n,b=*(cnt.points_begin());
    Vector_2 k2=::Vector_3to2(k,Vector_3(),x,y),b2=Vector_3to2(b);
    vec_type frac[2];
    Vector_2 surfp[2];
    int tl=base->TestLine(b2,k2,frac,surfp);
    if(tl<=1)
      return 0; // all outside, no intersection
    Plane_3 planes[3];
    planes[0].init(k, Vector_2to3(surfp[0]));
    planes[1].init(-k, Vector_2to3(surfp[1]));
    Polyhedron_3 poly(planes, planes+2, 0);
    return poly.TestContour(cnt, NULL, subcenter);
  }

  // контур проецируется на основание
  Contour_3to2<contour_t> cnt2(cnt,origin,x,y);
  Vector_2 subcenter2;
  vec_type area=base->TestContour(cnt2, NULL, subcenter ? &subcenter2 : NULL);
  area/=fabs(coef);
  if(area!=0 && subcenter){
    Vector_3 b=*(cnt.points_begin()); // начальная вершина контура
    vec_type hc=(b-origin)*n; // ее высота
    Vector_3 c=Vector_2to3(subcenter2)+hc*n; // переводит в 3D на уровне точки b
    Plane_3 pl(ncnt,b); // плоскость контура
    pl.TestRay(c,n,subcenter); // поднимаем в плоскость контура
  }
  return area;
}

template<class base_tt>
int Cylinder<base_tt>::TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
vector_t *surfp, vector_t *surfn, vec_type epsilon) const{
  Vector_3 v=dir%n;
//  if(acccomp(v.norm2()))
  if(acccomp(v.norm(),0.))
    return 0; // parallel to axis

  Vector_2 p2=Vector_3to2(p), dir2=::Vector_3to2(dir,0,x,y);
  Vector_2 surfn2[2];
  int tl=base->TestLine(p2,dir2,frac,NULL,surfn ? surfn2 : NULL);
  for(int i=0;i<tl;i++){
    if(surfp)
      surfp[i]=p+frac[i]*dir;
    if(surfn)
      surfn[i]=surfn2[i][0]*x+surfn2[i][1]*y;
  }
  return tl;
}

template<class base_tt>
int Cone<base_tt>::TestLine(const vector_t &p, const vector_t &dir, vec_type *frac,
vector_t *surfp, vector_t *surfn, vec_type epsilon) const{

  // we solve square equation for t:
  // (p + t*dir - origin, x)^2 + (p + t*dir - origin, y)^2 = (R * (1 - (p + t*dir - origin, n)/L)))^2

  return 1;
}

template<class reg_tt>
template<class contour_t>
vec_type StretchedRegion<reg_tt>::TestContour(const contour_t &cnt, VecContour<reg_tt::dimension> *subcont, typename reg_tt::vector_t *subcenter) const{
  VecContour<reg_tt::dimension> cnt1,subcont1;
  vector_t subcenter1;
  for(typename contour_t::point_iterator it=cnt.points_begin(),e=cnt.points_end();it!=e;++it)
    cnt1.add(basis_inv(*it-shift));
  vec_type area=reg->TestContour(cnt1,subcont ? &subcont1 : NULL,subcenter ? &subcenter1 : NULL);

  Vector_3 nn=cnt.GetNormVect(),x,y;
  set_perpendiculars(nn,x,y);
  Vector_3 x1=basis(x),y1=basis(y);
  vec_type change=sqrt((x%y).norm2()/(x1%y1).norm2());
  area*=change;

  if(subcont){
    subcont->Clear();
    for(int i=0;i<subcont1.GetNPoints();i++)
      subcont->add(shift+basis(subcont1.GetPoint(i)));
  }
  if(subcenter){
    *subcenter=shift+basis(subcenter1);
  }

  return area;
}

template<class reg_tt>
template<class contour_t>
vec_type ConfinedRegion<reg_tt>::TestContour(const contour_t &cnt, VecContour<> *subcont, Vector_3 *subcenter) const{
  VecContour<> subcnt;
  vec_type area=poly->TestContour(cnt, &subcnt, NULL);
  return area ? reg->TestContour(subcnt,subcont,subcenter) : 0;
}

template<class reg_tt>
int ConfinedRegion<reg_tt>::TestLine(const Vector_3 &p, const Vector_3 &dir, vec_type *frac,
Vector_3 *surfp, Vector_3 *surfn, vec_type epsilon) const{
  vec_type frac2[2][2];
  vector_t surfp2[2][2],surfn2[2][2];
  int res[2];
  bool in[2];
  Region<dimension> *r[2]={reg.ptr(),poly.ptr()};
  for(int i=0;i<2;i++){
    res[i]=r[i]->TestLine(p,dir,frac2[i],surfp2[i],surfn2[i]);
    in[i]=r[i]->TestPoint(p);
  }
  if(in[0] && in[1]){ // inside
    vec_type min=-VEC_INFTY, max=VEC_INFTY;
    int ii[2]={-1,-1,},ij[2]={-1,-1};
    for(int i=0;i<2;i++){
      for(int j=0;j<res[i];j++){
        if(frac2[i][j]>min && frac2[i][j]<0){
          min=frac2[i][j];
          ii[0]=i, ij[0]=j;
        }
        if(frac2[i][j]<max && frac2[i][j]>0){
          max=frac2[i][j];
          ii[1]=i, ij[1]=j;
        }
      }
    }
    int ind=0;
    for(int i=0;i<2;i++){
      if(ii[i]>=0){
        frac[ind]=frac2[ii[i]][ij[i]];
        if(surfp)
          surfp[ind]=surfp2[ii[i]][ij[i]];
        if(surfn)
          surfn[ind]=surfn2[ii[i]][ij[i]];
        ind++;
      }
    }
    return ind;
  }
  else if(!in[0] && !in[1]){ // outside both regions... not tested!
    for(int i=0;i<2;i++){
      if(res[i]<=0) // line does not intersect at least one region
        return 0;
    }
    if((frac2[0][0]>=0 && frac2[1][0]<=0) || (frac2[1][0]>=0 && frac2[0][0]<=0)) // regions are at different sides of line
      return 0;

    int open = fabs(frac2[1][0]) > fabs(frac2[0][0]); // far opening region
    vec_type open_val = fabs(frac2[open][0]);
    int close=-1;
    vec_type close_val=VEC_INFTY;
    for(int i=0;i<2;i++){
      if(res[i]==2){
        if(fabs(frac2[i][1]) < close_val){
          close_val=fabs(frac2[i][1]);
          close=i;
        }
      }
    }
    if(open_val>=close_val)
      return 0;
    frac[0]=frac2[open][0];
    if(surfp)
      surfp[0]=surfp2[open][0];
    if(surfn)
      surfn[0]=surfn2[open][0];

    if(close==-1)
      return 1;

    frac[1]=frac2[close][1];
    if(surfp)
      surfp[1]=surfp2[close][1];
    if(surfn)
      surfn[1]=surfn2[close][1];

    return 2;
  }
  // other cases still are not programmed
  return 0;
}
/*
template<class reg_tt>
vec_type ConfinedRegion<reg_tt>::TestEdge(const Vector_3 &p1, const Vector_3 &p2, Vector_3 *surfp, Vector_3 *surfn) const{
  Vector_3 surfp2[2], surfn2[2];
  vec_type res[2];
  res[0]=reg->TestEdge(p1, p2, surfp2, surfn2);
  res[1]=poly->TestEdge(p1, p2, surfp2+1, surfn2+1);
  int ires;
  if (res[0]==1 && res[1]==1)
    return 1;
  if (res[0]==-1 || res[1]==-1)
    return -1;
  if (res[0]==1)
    ires=1;
  else if (res[1]==1)
    ires=0;
  else {
    if (reg->TestPoint(surfp2[1]))
      ires=1;
    else
      ires=0;
  }
  if (surfp)
    *surfp=surfp2[ires];
  if (surfn)
    *surfn=surfn2[ires];
  return res[ires];
}
*/

template<class inp_it>
int MakeBoxDomains(const Box &box,int nproc, int nx, int ny, int nz, inp_it result){
  int ndir[3]={nx,ny,nz}, nspl[3]={nx,ny,nz};
  int np=(nproc>0 ? nproc : -nproc);
  int startdir=2;
  int dd, curdir=-1;
  // getting minimal splits
  int nauto=0, autod[3]={-1, -1, -1}, ndef=0;
  for(dd=0;dd<3;dd++){
    int dir=(dd+startdir)%3;
    if(ndir[dir]>0 && ndir[dir]>curdir)curdir=dir;
    else if(ndir[dir]<0){
      autod[nauto++]=dir;
    }
    if(ndir[dir]==0)ndef++;
  }
  if(curdir<0 && nauto>0){
    if(nauto==1){
      ndir[autod[0]]=np;
      curdir=autod[0];
    }
    else{  // auto-finding direction and number of slices to split
      int i;
      vec_type lprod=1.;
      for(i=0;i<nauto;i++){
        lprod*=box.GetSize()[autod[i]];
      }
      vec_type k=pow((vec_type)np/lprod,vec_type(1)/vec_type(nauto));
      int mnp=0, mdir=0;
      for(i=0;i<nauto;i++){
        int tnp=(int)(box.GetSize()[autod[i]]*k);
        if(mnp<tnp){
          mnp=tnp;
          mdir=i;
        }
      }
      if(mnp>0){
        ndir[autod[mdir]]=mnp;
        curdir=autod[mdir];
      }
    }
  }
  if(curdir>=0){ // splitting by minimal dir, auto split is ignored
    // number of domains
    int ndom=0;
    // determining sub-nps
    // THIS IS NP SWITCH : !!!
    if(np<ndir[curdir])ndir[curdir]=np;
    if(np>ndir[curdir] && ndef==2){ // the last axis
      ndir[curdir]=np;
    }

    int snp=(np>ndir[curdir] ? np : ndir[curdir]);
    int mnp=snp/ndir[curdir]; // minimal np
    int dntot=snp%ndir[curdir]; // difference
    vec_type c0=(vec_type)mnp/snp, c1=(vec_type)(mnp+1)/snp;
    Vector_3 p1=box.get_p1();
    Vector_3 p2=box.get_p2();
    
    vec_type l=box.GetSize()[curdir];
    int i;
    for(i=0;i<ndir[curdir];i++){
      vec_type dl=(dntot>0 ? l*c1 : l*c0);
      int cnp=(dntot>0 ? mnp+1 : mnp);
      if(cnp<=0)break;
      dntot--;
      p2[curdir]=p1[curdir]+dl;
      Box subdiv(p1,p2);
      nspl[curdir]=0; // this direction is processed
      int ddom=MakeBoxDomains(subdiv,(nproc>0 ? cnp: -cnp),nspl[0],nspl[1],nspl[2],result);
      if(nproc>0)result+=ddom;
      ndom+=ddom;
      p1[curdir]=p2[curdir];
    }
    return ndom;
  }
  /// called with  (0,0,0) split: put this box into the set
  if(nproc>0)
    *result=box;
  return 1;   
}
