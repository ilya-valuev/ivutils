# ifndef __PPART_H
# define __PPART_H

# include "partition.h"
# include "vector_3.h"

class ParticleDistribution: public SpaceDistribution {
protected:
  const Vector_3 *Coords;
  long nc;
  long *indarr;

  struct iGroup {
    void *pptr;
    long ng;
    long *indarr;
    int operator==(const iGroup &o){
      if(pptr==o.pptr && ng==o.ng &&  indarr==o.indarr)return 1;
      else return 0;
    }
  };

  lList<iGroup> constr;
  void init(long ncoord,const Vector_3 *Coordarr,long *inda);

  // used by split creator
  virtual SpaceDistribution* ExtractSub() const {
    SpaceDistribution sd(*this);
    sd.Refine();
    return new ParticleDistribution(indarr,nc,Coords,sd);
  }

public:

  ParticleDistribution(long ncoord,const Vector_3 *Coordarr){
    init(ncoord,Coordarr,NULL);
  }

  ParticleDistribution(long *inda,long ncoord,const Vector_3 *Coordarr,
                       const SpaceDistribution& sd):SpaceDistribution(sd){
    init(ncoord,Coordarr,inda);
  }



  virtual void *AddPlaneConstraint(int ic,vec_type val,int dir);

  virtual int RemovePlaneConstraint(void *ptr);

  virtual void* ReversePlaneConstraint(void *ptr);

  virtual void GetCoordDistr(int ic, Distribution *dens);

  virtual ~ParticleDistribution(){
    delete indarr;
  }

  virtual int write(char *filename){
    FILE *f=Err_fopen(filename,"w+t");
    long i;
    fprintf(f,"#Npart: %ld\n",nc);
    for(i=0;i<nc;i++){
      fprintf(f,"%g %g %g\n",Coords[indarr[i]][0],Coords[indarr[i]][1],Coords[indarr[i]][2]);
    }
    fclose(f);
    return 1;
  }


};




# endif