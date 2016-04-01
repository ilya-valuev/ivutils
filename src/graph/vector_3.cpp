/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.13 $
 *   $Date: 2012/06/29 10:50:13 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/graph/vector_3.cpp,v 1.13 2012/06/29 10:50:13 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/graph/vector_3.cpp,v $
$Revision: 1.13 $
$Author: valuev $
$Date: 2012/06/29 10:50:13 $
*/
/*s****************************************************************************
 * $Log: vector_3.cpp,v $
 * Revision 1.13  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.14  2012/04/13 12:41:24  valuev
 * moved stdlib.h back (NULL is defined there under llinux)
 *
 * Revision 1.11  2009/07/24 05:08:21  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.10  2009/06/01 13:01:52  valuev
 * Added ShearBox
 *
 * Revision 1.9  2009/01/30 13:54:24  valuev
 * restructured as a library
 *
 * Revision 1.9  2008/07/02 13:11:19  valuev
 * new C60+O2 experiments
 *
 * Revision 1.8  2008/04/15 13:11:16  valuev
 * Added antisymmetrized wave packets
 *
 * Revision 1.7  2008/01/29 23:01:20  valuev
 * Added VASP interface
 *
 * Revision 1.6  2007/07/09 21:27:49  valuev
 * plasma with wave packets
 *
 * Revision 1.6  2007/04/17 10:50:40  valuev
 * Added new contour path models: volume and tensor
 *
 * Revision 1.5  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.4  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
 * Revision 1.3  2006/10/27 20:41:02  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.2  2006/10/24 11:42:16  valuev
 * added emCorrectionSet class for packing corrections
 *
 * Revision 1.1  2006/08/24 12:26:16  valuev
 * added missing files
 *
 * Revision 1.2  2006/08/08 12:22:43  valuev
 * Added geometry
 *
 * Revision 1.1  2005/12/09 21:06:46  valuev
 * Added neighbour list to mdPotential interface.
 * Added missing files to ivutils directory.
 * Added mdtutorial and step1 project
 *
 *
*******************************************************************************/
//# include <iostream.h>
# include <stdio.h>


#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

# include "vector_3.h"


vec_type dist_av(Vector_3 *va1,Vector_3 *va2,int n){
 vec_type d;
 Vector_3 v;
 int i;

 d=0.;
 for(i=0;i<n;i++){
  v=va2[i]-va1[i];
  d+=v.norm();
 }
 d/=n;
 return d;
}

vec_type dist_max(Vector_3 *va1,Vector_3 *va2,int n){
 vec_type d,m;
 Vector_3 v;
 int i;

 m=0.;
 for(i=0;i<n;i++){
  v=va2[i]-va1[i];
  d=v.norm();
  if(m<d)m=d;
 }
 return m;
}


vec_type diff_av(Vector_3 *va1,Vector_3 *va2,int n, int *minind, int *maxind){
  int i;
  vec_type sum=0, nmax=-1, nmin=1e32, cnrm;
  int imin, imax;
  for(i=0;i<n;i++){
    if(va2)cnrm=(va1[i]-va2[i]).norm();
    else cnrm=va1[i].norm();
    if(cnrm>nmax){
      nmax=cnrm;
      imax=i;
    }
    if(cnrm<nmin){
      nmin=cnrm;
      imin=i;
    }
    sum+=cnrm;
  }
  sum/=n;
  if(minind)*minind=imin;
  if(maxind)*maxind=imax;
  return sum;
}


Vector_3 GetScope(const Vector_3 *varr,long n,Vector_3* cube1,Vector_3* cube2){
 Vector_3 center(varr[0]);
 *cube1=varr[0];
 *cube2=varr[0];
 for(long i=1;i<n;i++){
  center+=varr[i];
  for(int j=0;j<3;j++){
   if((*cube1)[j]>varr[i][j])(*cube1)[j]=varr[i][j];
   if((*cube2)[j]<varr[i][j])(*cube2)[j]=varr[i][j];
  }
 }
 return center/n;
}

Vector_3 GetIScope(const Vector_3 *varr,long *indarr,long n,Vector_3* cube1,Vector_3* cube2){
 Vector_3 center(varr[indarr[0]]);
 *cube1=varr[indarr[0]];
 *cube2=varr[indarr[0]];
 for(long i=1;i<n;i++){
  center+=varr[indarr[i]];
  for(int j=0;j<3;j++){
   if((*cube1)[j]>varr[indarr[i]][j])(*cube1)[j]=varr[indarr[i]][j];
   if((*cube2)[j]<varr[indarr[i]][j])(*cube2)[j]=varr[indarr[i]][j];
  }
 }
 return center/n;
}


Vector_3 GetIScopei(const Vector_3 *varr,int *indarr,int n,Vector_3* cube1=NULL,Vector_3* cube2=NULL){
 Vector_3 center(varr[indarr[0]]);
 if(cube1)*cube1=varr[indarr[0]];
 if(cube2)*cube2=varr[indarr[0]];
 for(long i=1;i<n;i++){
  center+=varr[indarr[i]];
  for(int j=0;j<3;j++){
    if(cube1){
      if((*cube1)[j]>varr[indarr[i]][j])(*cube1)[j]=varr[indarr[i]][j];
    }
    if(cube2){
      if((*cube2)[j]<varr[indarr[i]][j])(*cube2)[j]=varr[indarr[i]][j];
    }
  }
 }
 return center/n;
}

void clear_vecarri(int n,Vector_3 *vec, int *ind){
  if(!vec)return;
  int i, c;
  for(i=0;i<n;i++){
    if(ind)c=ind[i];
    else c=i;
    vec[c]=0.;
  }
}


Vector_3 FindPerp(const Vector_3 &vAB){
  Vector_3 res;
  // finding max. comp
  double vmax=fabs(vAB[0]);
  double vmin=vmax;
  int maxind=0, minind=0,i;
  for(i=1;i<3;i++){
    double v=fabs(vAB[i]);
    if(vmax<v){
      vmax=v;
      maxind=i;
    }
    if(vmin>v){
      vmin=v;
      minind=i;
    }
  }
  double s=0;
  if(minind==maxind){// zero vector
    res=Vector_3(1,1,1);
  }
  else{
    res[minind]=1.;
    res[maxind]=-vAB[minind]/vAB[maxind];
  }
  res.normalize();
  return res;
}

///\en Returns \a true if absolute value of \a a is at least \a factor times less than
///    the absolute value of \a b.
template<class T>
bool much_less(const T& a, const T& b, const T &factor){
  if(fabs(a)<fabs(b))
    return fabs(a*factor/b)<1 ? true : false;
  else
    return false;
}

Vector_3 Reflect(Vector_3& ini, double t, Vector_3 &vdir, double *box, int flag, const Vector_3 &force){
  Vector_3 dxforce=force*(t*t/2);
  Vector_3 vres, vini=ini;
  double tc, tcross;
  int icross;
  int i, nc;
  do{
    Vector_3 dxforce=force*(t*t/2);
    vres=ini+vdir*t+dxforce;
    // counting crossings
    nc=0;
    tcross=-1;
    for(i=0;i<3;i++){ // with 0-side
      if( (flag&(0x1<<i)) && vres[i]<0){
        double dx=0-vini[i];
        if(much_less(force[i]*dx,vdir[i]*vdir[i],1e5)) // linear equation
          tc=dx/vdir[i];
        else { // quadratic
          if(much_less(dxforce[i],dx,1e5)) // force is too small
            continue;
          double det=vdir[i]*vdir[i]+2*force[i]*dx;
          if(det<0)
            continue;
          tc=(sqrt(det)-vdir[i])/force[i];
        }
        if(tcross<0 || tc<tcross){
          icross=i;
          tcross=tc;
        }
        nc++;
      }
    }
    for(i=0;i<3;i++){ // with box-side
      if( (flag&(0x1<<i)) && vres[i]>box[i]){
        double dx=box[i]-vini[i];
        if(much_less(force[i]*dx,vdir[i]*vdir[i],1e5)) // linear equation
          tc=dx/vdir[i];
        else { // quadratic
          if(much_less(dxforce[i],dx,1e5)) // force is too small
            continue;
          double det=vdir[i]*vdir[i]+2*force[i]*dx;
          if(det<0)
            continue;
          tc=(sqrt(det)-vdir[i])/force[i];
        }
        //tc=(box[i]-vini[i])/vdir[i];
        if(tcross<0 || tc<tcross){
          icross=i;
          tcross=tc;
        }
        nc++;
      }
    }
    // reflecting
    if(tcross>0 && tcross<t){
      vini=vini+vdir*tcross+force*(tcross*tcross/2);
      t-=tcross;
      vdir[icross]=-vdir[icross];
      continue;
    }
  }while(0);
  return vres;
}

Vector_3 randdir(){
  vec_type xi1 = 2.*((vec_type)rand()) / RAND_MAX-1.;
  vec_type xi2 =   ((vec_type)rand()) / RAND_MAX;
  vec_type r1  = sqrt (1.-xi1*xi1);
  return  Vector_3(r1*cos(2.*M_PI*xi2),r1*sin(2.*M_PI*xi2),xi1);
}






/*
void main(){
 Vector_3 vector,v1(0,34,6),v2(1,0,0);
 float c;
 vector.print();
 v1.print();
 vector=v1;
 vector.print();
 vector=v1+v2;
 vector.print();
 v1.print();

 c=v1*vector;
 cout<<c<<"\n";
 v1*=c;
 v1.print();
 c=v1.norm();
 v1.normalize();
 v1=v1%v1;
 v1.print();
 _asm{
  mov ax,10
  mov dx,ax
 }
} */







