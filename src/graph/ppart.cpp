# include "ppart.h"


void ParticleDistribution::init(long ncoord,const Vector_3 *Coordarr,long *inda){
  indarr= new long[ncoord];
  if(!indarr)fatal_error("ParticleDistribution: MAE!\n");
  long i;
  if(!inda){
    for(i=0;i<ncoord;i++)indarr[i]=i;
  }
  else{
    for(i=0;i<ncoord;i++)indarr[i]=inda[i];
  }

  Vector_3 pmin, pmax;
  GetIScope(Coordarr,indarr,ncoord,&pmin,&pmax);
  SpaceRegion lim(pmin,pmax);
  SetLimitRegion(&lim);
  Coords=Coordarr;
  nc=ncoord;
}

void *ParticleDistribution::AddPlaneConstraint(int ic,vec_type val,int dir){
  void *ptr=AddCoordPlane(ic,val,dir);
  // sorting array
  long i, top=nc;
  if(dir>0){  // leaving all less
    for(i=0;i<top;i++){
      if(Coords[indarr[i]][ic]>=val){
        // swaping
        int tmp=indarr[top-1];
        indarr[top-1]=indarr[i];
        indarr[i]=tmp;
        top--;
        i--;
      }
    }
  }
  else{
    for(i=0;i<top;i++){  // leaving all greater or equal
      if(Coords[indarr[i]][ic]<val){
        // swaping
        int tmp=indarr[top-1];
        indarr[top-1]=indarr[i];
        indarr[i]=tmp;
        top--;
        i--;
      }
    }
  }
  iGroup g;
  g.pptr=ptr;
  g.ng=nc;
  g.indarr=indarr;
  constr.Append(g);
  nc=top;
  return ptr;
}

void* ParticleDistribution::ReversePlaneConstraint(void *ptr){
  constr.Rewind();
  while(!constr.Done()){
    if(constr.GetCur().pptr==ptr){
      iGroup g=constr.GetCur();

      // removing till the end
      long rmi=constr.GetCurPos();
      constr.Next();
      while(!constr.Done()){ // removing all added afterwards
        RemoveCoordPlane(constr.GetCur().pptr);
        constr.Next();
      }
      constr.Remove(rmi,-1); // removing entries added afterwards and current

      g.pptr=SpaceDistribution::ReversePlaneConstraint(ptr);

      if(indarr==g.indarr){ // base was not changed
        indarr=g.indarr+nc; //changing base
      }
      else{  // base was already changed
        indarr=g.indarr; // returning base
      }
      nc=g.ng-nc; // changing number
      constr.Append(g);

      return g.pptr;
    }
    constr.Next();
  }
  return NULL;
}

int ParticleDistribution::RemovePlaneConstraint(void *ptr){
  constr.Rewind();
  while(!constr.Done()){
    if(constr.GetCur().pptr==ptr){
      iGroup g=constr.GetCur();
      // removing till the end
      long rmi=constr.GetCurPos();

      do{ // removing all added afterwards and current
        RemoveCoordPlane(constr.GetCur().pptr);
        constr.Next();
      }while(!constr.Done());
      constr.Remove(rmi,-1); // removing entries added afterwards and current

      nc=g.ng;   // restoring variables
      indarr=g.indarr;
      return 1;
    }
    constr.Next();
  }
  return 0;
}



void ParticleDistribution::GetCoordDistr(int ic, Distribution *dens){
  dens->clear();
  long i;
  for(i=0;i<nc;i++){
    if(Test((Vector_3 *)&Coords[indarr[i]])>=0){
      dens->point(Coords[indarr[i]][ic],1);
    }
  }
}