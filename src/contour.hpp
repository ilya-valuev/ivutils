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
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/contour.hpp,v 1.3 2012/06/29 10:50:13 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/contour.hpp,v $
$Revision: 1.3 $
$Author: valuev $
$Date: 2012/06/29 10:50:13 $
*/
/*s****************************************************************************
 * $Log: contour.hpp,v $
 * Revision 1.3  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.16  2011/05/12 21:35:32  lesha
 * small_coord is changed
 *
 * Revision 1.15  2011/01/27 18:14:54  lesha
 * small_coord value is changed in ProjectSimplex
 *
 * Revision 1.14  2010/06/08 22:35:04  lesha
 * *** empty log message ***
 *
 * Revision 1.13  2010/06/08 12:09:46  lesha
 * *** empty log message ***
 *
 * Revision 1.12  2010/01/18 12:46:41  lesha
 * *** empty log message ***
 *
 * Revision 1.11  2010/01/18 11:47:13  lesha
 * *** empty log message ***
 *
 * Revision 1.10  2010/01/03 13:51:52  lesha
 * TestLine is added
 *
 * Revision 1.9  2010/01/02 22:22:58  lesha
 * StretchedRegion is added, Contour::GetCenter if fixed, 2DDumping is added etc.
 *
 * Revision 1.8  2010/01/01 20:40:24  lesha
 * Region, Box_N, Sphere_N, Basis_N, MonteCarloTestContour are added
 *
 * Revision 1.7  2009/12/01 22:42:54  lesha
 * Area() is moved into Contour. aggregate_it is moved to utiltl.h
 *
 * Revision 1.6  2009/12/01 18:46:07  lesha
 * contour is modified
 *
 * Revision 1.5  2009/09/05 14:47:53  lesha
 * DumpPolyhedron and GetFaces are modified
 *
 * Revision 1.4  2009/06/21 10:38:47  lesha
 * region_2.h is added
 *
 * Revision 1.3  2009/06/12 20:17:12  lesha
 * ProjectSimplex is fixed
 *
 * Revision 1.2  2009/05/22 08:53:52  valuev
 * new uiExperiment interface
 *
 * Revision 1.1  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.8  2008/01/23 18:14:17  lesha
 * template class instantiation is corrected
 *
 * Revision 1.7  2008/01/23 17:10:49  lesha
 * ...
 *
 * Revision 1.6  2008/01/23 17:08:12  lesha
 * make gcc compilable
 *
 * Revision 1.5  2008/01/22 10:16:27  lesha
 * unnecessary include files are removed
 *
 * Revision 1.4  2008/01/06 03:25:06  lesha
 * Implementation part is moved from *.h to *.cpp
 *
 * Revision 1.3  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/

/// \file \brief Template function definitions for contour.h

#include "contour.h"
#include "pencil.h"
#include "utiltl.h"

extern pencil<Vector_3> _PStmpv;
extern pencil<int> _PStmps;

template<class point_it,int N>
Vector_Nt<vec_type,N> GetBoundingBox(point_it it, point_it end,Vector_Nt<vec_type,N>* cube1,Vector_Nt<vec_type,N>* cube2){
  Vector_Nt<vec_type,N> center(*it);
  *cube1=*it;
  *cube2=*it;
//  it++;
  int n=0;
  for(;it!=end;it++){
    Vector_Nt<vec_type,N> v=*it;
    center+=v;
    for(int j=0;j<N;j++){
      if((*cube1)[j]>v[j])(*cube1)[j]=v[j];
      if((*cube2)[j]<v[j])(*cube2)[j]=v[j];
    }
    n++;
  }
  if(n)center/=n;
  return center;
}

template<class point_it,class edge_it,class edge_t>
vec_type Contour<point_it,edge_it,edge_t>::Area() const{
  point_iterator it=points_begin(), e=points_end();
  if(!(it!=e))return 0.; // 0 points in the contour?
  vector_t p[3];
  p[0]=*it++;
  if(!(it!=e))return 0.; // 1 point in the contour?
  p[1]=*it++;
  if(!(it!=e))return 0.; // 2 points in the contour?
  p[2]=*it;
  vector_t v1=p[0]-p[1];
  vector_t v2=p[2]-p[1];

  vec_type area=0;
  do{
    area+=vec_area(v1, v2)/2;
    ++it;
    if(!(it!=e))break;
    v1=p[0]-p[2];
    v2=*it-p[2];
    p[2]=*it;
  }while(1);

  return area;
}

template<class point_it,class edge_it,class edge_t>
typename point_it::value_type Contour<point_it,edge_it,edge_t>::GetCenter() const{
  point_iterator it=points_begin(),e=points_end();
  if(!(it!=e))return 0; // 0 points in the contour?
  vector_t c=*it;
  ++it;
  if(!(it!=e))return 0; // 1 point in the contour?
  vector_t v1=*it-c;
  ++it;
  if(!(it!=e))return 0; // 2 points in the contour?
  vector_t v2=*it-c;
  ++it;

  vec_type area=vec_area(v1,v2)/2;
  vector_t center=area*(v1+v2)/3; // medians intersection

  while(it!=e){
    v1=v2;
    v2=*it-c;
    vec_type darea=vec_area(v1,v2)/2;
    area+=darea;
    center+=darea*(v1+v2)/3;
    ++it;
  }
  return c+center/area;
}

template<class point_it,class edge_it,class edge_t>
bool Contour_N<point_it,Vector_2,edge_it,edge_t>::TestPoint(const Vector_2 p) const{
  point_it it=points_begin(), e=points_end();
  Vector_2 a=*it-p,b;
  int rotate=0;
  bool first=true;
  while(first){
    ++it;
    if(!(it!=e)){
      it=points_begin();
      first=false;
    }
    b=*it-p;
    vec_type cur=a[0]*b[1]-a[1]*b[0];
    if(cur){
      if(!rotate)
        rotate=cur>0 ? 1 : -1;
      else if (rotate!=(cur>0 ? 1 : -1)) 
        return false;
    }
    a=b;
  }
  return true;
}

template<class point_it,class edge_it,class edge_t>
bool Contour_N<point_it,Vector_3,edge_it,edge_t>::TestEdge(const Vector_3 &v1, const Vector_3 &v2, Vector_3 *cross) const{
  Plane_3 pl=GetPlane();
  Vector_3 cr;
  if(!pl.TestEdge(v1,v2,&cr))
    return false;
  point_it it=points_begin(), e=points_end();
  edge_it ite=base_t::edges_begin(), ee=base_t::edges_end();
  Vector_3 nv=GetNormVect();
  for(; it!=e; ++it,++ite){
    Vector_3 r=cr-*it;
    if( ((*ite).rel()%r)*nv <0)
      return false; //outside 
  }
  if(cross)
    *cross=cr;
  return true;
}

template<class point_it,class edge_it,class edge_t>
Vector_3 Contour_N<point_it,Vector_3,edge_it,edge_t>::GetNormVect() const {
  point_it it=points_begin(), e=points_end();
  Vector_3 p[3];
  int i=0;
  for(; it!=e && i<3; it++, i++)
    p[i]=*it;
  if(i<3)return Vector_3(0,0,0);
  Vector_3 v1=p[0]-p[1], v2=p[2]-p[1];
  return v2%v1;
}

/*
template<class contour_t>
int project_contour3(contour_t &cont3, VecContour<2> &cont2, Vector_3 c, Vector_3 x, Vector_3 y){
  cont2.Clear();
  for(typename contour_t::point_iterator it=cnt.points_begin(), e=cnt.points_end();it!=e;++it){
    Vector_3 p=*it;
    Vector_2 p2=Vector_3to2(p,c,x,y);
    cont2.add(p2);
  }
  return 1;
}
*/
struct sedge_t{
  Vector_3 p[2];
  int k[2];
  int i;
  sedge_t(int kl=0, int kr=0, int k0=0, Vector_3 pl=0, Vector_3 pr=0): i(k0) {
    k[0]=kl;
    k[1]=kr;
    p[0]=pl;
    p[1]=pr;
  }
};

template<class plane_it, class contour_t>
int ProjectSimplex(const Plane_3& plane0, plane_it beg, plane_it end, contour_t &cnt, int nplanes){
  double small=VEC_ZERO, small_coord=1e-10, small_coord2=small_coord*small_coord;
  cnt.Clear();
  int neq=0;
  // counting
  
  if(nplanes<=0){
    for(plane_it it=beg;it!=end;){
      neq++;
      it++;
    }
  }
  else neq=nplanes;
  vector<sedge_t> edges;

  // lines
  _PStmpv.resize(neq);
  _PStmps.resize(neq);

  Vector_3 *peq= _PStmpv.get_ptr();
  int *state= _PStmps.get_ptr();
  Vector_3 n0;
  vec_type d0;
  plane0.get_coeff(n0,d0);
  // determinig which coord to substitute
  int im, io[2];
  vec_type vm=n0.maxabscoord(&im);
  io[0]=(im+1)%3;
  io[1]=(im+2)%3;
  // Делим на vm, после этого уравнение плоскости n0*point+d0=0 остается в силе
  // После деления получаем n0[im]=1
  n0/=vm;
  d0/=vm;
  // substituting
  int i=0;
  for(plane_it it=beg;it!=end;it++){
    Vector_3 ni;
    vec_type di;
    (*it).get_coeff(ni,di);
    // Найдем прямую, являющуюся пересечением двух плоскостей, заданных уравнениями:
    // n0[io[0]]*r[io[0]]+n0[io[1]]*r[io[1]]+r[im]+d0=0
    // ni[io[0]]*r[io[0]]+ni[io[1]]*r[io[1]]+ni[im]*r[im]+di>=0
    // Для этого домножим первое уравнение на ni[im] и подставим в первое уравнение r[im], выраженное из второго уравнения
    // Получим, что искомая прямая задается уравнением peq[i][0]*r2[0]+peq[i][1]*r2[1]+peq[i][2]>=0,
    // где r2[i]=r[io[i]], а r[im] ищется из уравнения плоскости (любой из двух)
    peq[i][0]=ni[io[0]]-ni[im]*n0[io[0]];
    peq[i][1]=ni[io[1]]-ni[im]*n0[io[1]];
    peq[i][2]=di-ni[im]*d0;
    state[i]=0;
    i++;
  }
  // solving simplex
  int res0=1, basic=0;
  int kstart=-1;

  for(i=0;i<neq;i++){
    // determining max coord
    int jm = (fabs(peq[i][0])>fabs(peq[i][1])) ? 0 : 1;
    int jo=(jm+1)%2; // complimentary
    if(fabs(peq[i][jm])<small){ // zero coefficients
      // плоскость plane0 лежит в полупространстве, выделяемом гранью многогранника
      if(peq[i][2]>=-small){ // >=0
        state[i]=-1; // skip
        continue; // trivial inequality
      }
      // наоборот, следовательно решения нет
      else{
        state[i]=-2; // contradictive
        res0=-1; // no solution
        break;
      }
    }
    vec_type tleft=-VEC_INFTY, tright=VEC_INFTY;
    int kright=-1;
    int kleft=-1;
    int res1=0;
    for(int k=0;k<neq;k++){ // comparing with others
      if(k==i || state[k]<0) // та же плоскость или плоскость, параллельная плоскости plane0
        continue;
      // Ищем пересечение двух прямых заданых уравнениями peq[j][jo]*r[jo]+peq[j][jm]*r[jm]+peq[j][2]>=0, j=i,k
      // Выражаем в первом уравнении r[jm]=-(peq[i][jo]*r[jo]+peq[i][2])/peq[i][jm]
      // и подставляем во второе уравнение:
      // (peq[k][jo]-peq[i][jo]*peq[k][jm]/peq[i][jm])*r[jo]+peq[k][2]-peq[i][2]*peq[k][jm]/peq[i][jm]>=0
      vec_type mult=peq[k][jm]/peq[i][jm];
      vec_type b=peq[k][2]-mult*peq[i][2]; // at+b>0
      vec_type a=peq[k][jo]-mult*peq[i][jo];
      // r[jo] принимаем за t
      // тогда a*t+b>=0
      if(fabs(a)<small){ // нет пересечения
        if(b>-small){
          continue; 
        }
        else{ // the current is trivial: move to kth
          int sign=Vector_2(&(peq[k][0]))*Vector_2(&(peq[i][0]))>0; // в одну ли сторону направлены нормали к прямым на плоскости
          if(!sign){ // полуплоскости, определяемые прямыми, не перекрываются
            state[i]=-2;
            res0=-1;
            break; // no solution
          }
          // else  just move to another one
          state[i]=-1; // прямая i лежит вне полуплоскости, выделяемой прямой k. Выбрасываем ее
          res1=-1; // i-th line is not in simplex
//          break;
        }
      }
      if(res1==-1)
        continue;
      vec_type t=-b/a;
      vec_type dt=t-tleft;
      // попутная прямая, двигаю t вправо
      if(a>0 && dt>-small){
        if(fabs(dt)>small){ 
          tleft=t;
          kleft=k;
        }
      }
      dt=t-tright;
      // встречная прямая, двигаю t влево
      if(a<0 && dt<small){
        if(fabs(dt)>small){
          tright=t;
          kright=k;
        }
      }
      dt=tleft-tright;
      if(dt>small){
        res1=-1; // на этой прямой нет отрезка
        break;
      }   
    }
    if(res0==-1)
      break;
    if(res1<0) // no edges found
      continue; 

    vec_type dt=tleft-tright;
    if(dt>-small) // tleft and tright are indistinguishable
      continue;
    
    if(kleft==-1 || kright==-1){ // no limits, solution is not a contour
      res0=2;
      break;
    }
    // OK, here tleft< tright
    // filling the points
    Vector_3 pointl, pointr;
    // Уравнение i-ой прямой на плоскости имеет вид peq[i][0]*point[0]+peq[i][1]*point[1]+peq[i][2]=0
    pointl[io[jo]]=tleft;
    pointl[io[jm]]=(-peq[i][2]-peq[i][jo]*tleft)/peq[i][jm];
    // Находим оставшуюся компоненту из уравнения плоскости n0*point+d0=0
    pointl[im]=-d0-n0[io[jo]]*pointl[io[jo]]-n0[io[jm]]*pointl[io[jm]];
    // Действуем аналогично
    pointr[io[jo]]=tright;
    pointr[io[jm]]=(-peq[i][2]-peq[i][jo]*tright)/peq[i][jm];
    pointr[im]=-d0-n0[io[jo]]*pointr[io[jo]]-n0[io[jm]]*pointr[io[jm]];

    // checking edge size
    if((pointl-pointr).norm2()<small_coord2)
      continue; 

    // inserting the edge
    edges.push_back(sedge_t(kleft,kright,i,pointl,pointr));

  }  // end plane loop

  if(res0!=1) // no contour
    return res0;
 
  // now forming the connections
  int icur=0, j0=1;
  vec_type drmin=0;
  int jmin=0;
  for(vector<sedge_t>::iterator it1=edges.begin();it1!=edges.end();){
    
    if(drmin<small_coord2) // direct connection
      cnt.SetPoint(icur++,it1->p[1-jmin]);
    else{ // indirect connection
      cnt.SetPoint(icur++,it1->p[jmin]);
      cnt.SetPoint(icur++,it1->p[1-jmin]);
    }
    it1->i=-1; // marking as done
    // finding the closest point
    jmin=0;
    drmin=VEC_INFTY;
    vector<sedge_t>::iterator itmin=edges.end();
    for(vector<sedge_t>::iterator it2=edges.begin();it2!=edges.end();++it2){
      if(it2->i<0)
        continue;
      for(int j=0;j<2;j++){
        vec_type dr=(it1->p[j0]-it2->p[j]).norm2();
        if(drmin>=dr){
          // checking for a clone
          if(dr<small_coord2 && (it1->p[1-j0]-it2->p[1-j]).norm2()<small_coord2){
            // the same edge
            it2->i=-1;
            break;
          }
          drmin=dr;
          itmin=it2;
          jmin=j;
        }
      }
    }
    it1=itmin;
    j0=1-jmin;
  }
  if(icur<3){
    cnt.Clear();
    res0=-3;
  }
  //else{ // adding the end point
  //  cnt.SetPoint(icur++,edges[0].p[0]);
  //}

  return res0;
}

template<class out_cont_t, class plane_it>
int GetFaces(out_cont_t &cont,plane_it beg,plane_it end){
  plane_it it=beg;
  typedef typename out_cont_t::value_type cont_t;
  typedef remove_if_it<plane_it> selector_it;
  
  int i=0;
  for(;it!=end;++it){
    plane_it sel=it;
    plane_it seln=sel;
    seln++;
    selector_it selb(beg,end,sel,seln);
    selector_it sele(end,end,end,end);
    cont_t cnt;
    ProjectSimplex(*it,selb,sele,cnt); // cut considered plane it from polyhedron
    //ProjectSimplex(*it,beg,end,cnt);
    if(cnt.GetNPoints()>0){
      cont.push_back(cnt);
      i++;
    }
  }
  return i;
}
