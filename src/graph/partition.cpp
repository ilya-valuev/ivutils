# include "partition.h"
//# include "vector_3.h"


// ---------- generic  box SpaceRegion --------------------------------

SpaceRegion::SpaceRegion(Vector_3 Point1,Vector_3 Point2){
  int i;
  for(i=0;i<3;i++){
   Point_min[i]=fmin(Point1[i],Point2[i]);
   Point_max[i]=fmax(Point1[i],Point2[i]);
  }
}

// tells whether the point is inside
// returns 1 if inside, 0 -- on the border, -1 --outside
int SpaceRegion::Test(Vector_3* Point) const {
  int i, bord=0;
  vec_type ret;
  for(i=0;i<3;i++){
    ret=(*Point)[i]-Point_min[i];
    if(ret<0)return -1;
    if(ret==0)bord++;
    ret=Point_max[i]-(*Point)[i];
    if(ret<0)return -1;
    if(ret==0)bord++;
  }
  if(bord>0)return 0;
  else return 1;
}


// tells whether the whole region is inside
// returns 1 if inside, 0 -- if crossing, -1 --outside
// attention: -2  if the decision is not implemented!
// generic function works with the enclosure box only!
int SpaceRegion::RegionTest(const SpaceRegion& region) const{
  Vector_3 v1=region.GetMin();
  int ret1=Test(&v1);
  v1=region.GetMax();
  int ret2=Test(&v1);

  if(ret1>=0 && ret2>=0) return 1; // all inside
  if(ret1==-1 && ret2 ==-1) return -1; // all outside
  return 0;
}

// ---------- generic  box SpaceRegion --------------------------------

//---------------CoordHalfspace----------------------------------------

CoordHalfspace::CoordHalfspace(int coordn,vec_type coordval,int sign):SpaceRegion(){
  if(coordn<0 || coordn>=3)fatal_error("CoordHalfspace: invalid coordinate number %d!\n",coordn);
  if(sign!=-1 && sign!=1)fatal_error("CoordHalfspace: invalid plane direction %d!\n",sign);
  ic=coordn;
  dir=sign;
  val=coordval;

  if(dir>0)Point_max[ic]=val;
  else Point_min[ic]=val;
}

//---------------CoordHalfspace----------------------------------------



//---------------------RegionSet-----------------------------------


Vector_3 *RegionSet::tstpoint=NULL;


// a is not yet added
void RegionSet::test_minmax(const SpaceRegion *a, long rgbefore){
  if(rgbefore<0)rgbefore=NRegs();
  int i;
  const Vector_3 &pmin=a->GetMin(), &pmax=a->GetMax();
  if(rgbefore==0){ // the set is yet void
    Point_min=(Vector_3)pmin;
    Point_max=(Vector_3)pmax;
  }
  else{
    for(i=0;i<3;i++){
      if(settype==SET_OR){
        if(Point_min[i]>pmin[i])Point_min[i]=pmin[i];
        if(Point_max[i]<pmax[i])Point_max[i]=pmax[i];
      }
      if(settype==SET_AND){
        if(Point_min[i]<pmin[i])Point_min[i]=pmin[i];
        if(Point_max[i]>pmax[i])Point_max[i]=pmax[i];
      }
    }
  }
}


void RegionSet::test_allminmax(){
  int before=0;
  if(NPlanes()>0){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++){
      CoordPlane[i][j].Rewind();
      while(!CoordPlane[i][j].Done()){
        test_minmax(CoordPlane[i][j].GetCur(),before);
        before=1;
        CoordPlane[i][j].Next();
      }
    }
  }
  OtherRegs.Rewind();
  while(!OtherRegs.Done()){
    test_minmax(OtherRegs.GetCur(),before);
    before=1;
    OtherRegs.Next();
  }
}

// PLANE SET MUST BE NON-VOID!
int RegionSet::planeOR(int *reso) const{
  int i;
  long ncr, ncl;

  if(*reso==-2)*reso=1;
  else if(*reso==0)return 0; // nothing will change

  vec_type minl,maxr;

  for(i=0;i<3;i++){      // leaving only [-INF, maxr] [minl,INF]
    ncr=CoordPlane[i][1].Count();
    ncl=CoordPlane[i][0].Count();
    if(ncr==0 && ncl==0)continue;


    if(ncr>0)maxr=CoordPlane[i][1].Get(ncr-1)->GetVal();
    else maxr=VEC_INFTY;
    if(ncl>0)minl=CoordPlane[i][0].Get(0)->GetVal();
    else minl=-VEC_INFTY;

    if((*tstpoint)[i]>maxr && (*tstpoint)[i]<minl){
      *reso=-1;
      return 1;
    }
    if((*tstpoint)[i]==maxr || (*tstpoint)[i]==minl )*reso=0;
  }
  return 0;
}

// PLANE SET MUST BE NON-VOID!
int RegionSet::planeAND(int *reso) const{
  int i;
  long ncr, ncl;

  if(*reso==-2)*reso=1;
  else if(*reso==-1)return 0; // definitely outside

  vec_type minr,maxl;
  for(i=0;i<3;i++){      // leaving only [maxl, minr]
    ncr=CoordPlane[i][1].Count();
    ncl=CoordPlane[i][0].Count();
    if(ncl==0 && ncr==0)continue;

    if(ncr>0)minr=CoordPlane[i][1].Get(0)->GetVal();
    else minr=VEC_INFTY;

    if(ncl>0)maxl=CoordPlane[i][0].Get(ncl-1)->GetVal();
    else maxl=-VEC_INFTY;

    if((*tstpoint)[i]<maxl || (*tstpoint)[i]>minr){
      *reso=-1;
      return 0;
    }
    if((*tstpoint)[i]==minr || (*tstpoint)[i]==maxl)*reso=0; // border

  }
  return 1;
}


void RegionSet::Refine(){
  int i;
  long j, ncr, ncl;
  if(settype==SET_OR){
    vec_type minl,maxr;
    for(i=0;i<3;i++){      // leaving only (maxr, minl)
      ncr=CoordPlane[i][1].Count();
      ncl=CoordPlane[i][0].Count();
      if(ncr==0 && ncl==0)continue;

      if(ncr>1)for(j=0;j<ncr-1;j++)CoordPlane[i][1].Remove(0);
      if(ncr>0)maxr=CoordPlane[i][1].Get(0)->GetVal();
      else maxr=VEC_INFTY;

      //for(j=1;j<ncl;j++)CoordPlane[i][0].Remove(j);
      if(ncl>1)CoordPlane[i][0].Remove(1,-1);
      if(ncl>0)minl=CoordPlane[i][0].Get(0)->GetVal();
      else minl=-VEC_INFTY;

      if(ncr>0)Point_min[i]=-VEC_INFTY;
      else Point_min[i]=minl;
      if(ncl>0)Point_max[i]=VEC_INFTY;
      else Point_max=maxr;
    }
  }
  if(settype==SET_AND){
    vec_type minr,maxl;
    for(i=0;i<3;i++){      // leaving only (maxl, minr)
      ncr=CoordPlane[i][1].Count();
      ncl=CoordPlane[i][0].Count();
      if(ncl==0 && ncr==0)continue;

      //for(j=1;j<ncr;j++)CoordPlane[i][1].Remove(j);
      if(ncr>1)CoordPlane[i][1].Remove(1,-1);
      if(ncr>0)minr=CoordPlane[i][1].Get(0)->GetVal();
      else minr=-VEC_INFTY;

      for(j=0;j<ncl-1;j++)CoordPlane[i][0].Remove(0);
      if(ncl>0)maxl=CoordPlane[i][0].Get(0)->GetVal();
      else maxl=VEC_INFTY;

      Point_min[i]=maxl;  // may be void
      Point_max[i]=minr;
    }
  }
  nplanes=0;
  for(i=0;i<3;i++)for(j=0;j<2;j++)nplanes+=CoordPlane[i][j].Count();
  test_allminmax();
}

//---------------------RegionSet-----------------------------------



//---------------------SpaceDistribution-----------------------------------

vec_type SpaceDistribution::locate_midpoint(TableFunction *dens,int n1, int n2){
  float integ=dens->integral();
  integ*=n1;
  integ/=n1+n2;
  return dens->integral_reaches(integ);
}

int SpaceDistribution::Subdivide(int nparts,vec_type avlength, RegionSet *Result){
  if(nparts<=0)fatal_error("SpaceDistribution.Subdivide: invalid number of parts:%d!\n",nparts);
  if(nparts==1){ // generate final set
    SpaceRegion *Set=ExtractSub();
    Result->AddRegion(Set,1);
    return 1;
  }

  int nsub1=nparts/2;
  int nsub2=nsub1+nparts%2;

  Vector_3 vmin, vmax;
  vmin=GetFuncMin();
  vmax=GetFuncMax();


  float cstep=avlength/10;
  int i;

  int imid, sw=0;
  float xt,xmid=0, mindens=-1,densv;

  for(i=0;i<3;i++){
    int nfsteps= (vmax[i]-vmin[i]+2*avlength)/cstep+2;
    if(nfsteps<=0)continue;
    Distribution cdistr(vmin[i]-avlength,vmax[i]+avlength,nfsteps);
    GetCoordDistr(i,&cdistr);
    Smooth(cdistr.get_func(),10);
    xt=locate_midpoint(cdistr.get_func(),nsub1,nsub2);
    densv=cdistr(xt);
    if(mindens<0 || densv<mindens){
      imid=i;
      xmid=xt;
      mindens=densv;
      sw=0;
    }
    if(nsub2!=nsub1){ // changing direction
      xt=locate_midpoint(cdistr.get_func(),nsub2,nsub1);
      densv=cdistr(xt);
      if(mindens<0 || densv<mindens){
        imid=i;
        xmid=xt;
        mindens=densv;
        sw=1;
      }
    }
  }
  if(mindens<0){
    msg_error("SpaceDistribution.Subdivide: void region?\n");
    return 0;
  }

  void *ind=AddPlaneConstraint(imid,xmid,1);
  if(!Subdivide((sw ? nsub2:nsub1),avlength,Result))return 0;
  //RemovePlaneConstraint(ind);
  ind=ReversePlaneConstraint(ind);

  //ind=AddPlaneConstraint(imid,xmid,-1);
  if(!Subdivide((sw ? nsub1:nsub2),avlength,Result))return 0;
  RemovePlaneConstraint(ind);

  return 1;
}

//---------------------SpaceDistribution-----------------------------------
