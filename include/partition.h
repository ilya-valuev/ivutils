# ifndef PARTITION_H
# define PARTITION_H


# include "common.h"
# include "vector_3.h"
# include "listiter.h"
# include "statist.h"

// this is the first start of MPI program for modelling Si coulomb
// explosion with STW potential   03.06.2002


// generic class defining the region in 3D space
class SpaceRegion{
protected:
   // the space region limiting the region
  Vector_3 Point_min, Point_max;
  //int type;
  //int subtype;

public:
  // constructor for the whole space (default)
  SpaceRegion(){
    Point_min=Vector_3(-VEC_INFTY,-VEC_INFTY,-VEC_INFTY);
    Point_max=Vector_3(VEC_INFTY,VEC_INFTY,VEC_INFTY);
  }

  // default constructor for a box
  SpaceRegion(Vector_3 Point1,Vector_3 Point2);

  // copy creation
  SpaceRegion(const SpaceRegion& other):Point_min(other.Point_min),Point_max(other.Point_max){}


  // clone function
  virtual SpaceRegion *Clone() const {
    return new SpaceRegion(*this);
  }

   // getting (box!) limits
  Vector_3 GetMin() const { return Point_min; };
  Vector_3 GetMax() const { return Point_max; };

  // tells whether the point is inside
  // returns 1 if inside, 0 -- on the border, -1 --outside
  virtual int Test(Vector_3* Point) const;

  // tells whether the whole region is inside
  // returns 1 if inside, 0 -- if crossing, -1 --outside
  // attention: -2  if the decision is not implemented!
  // generic function works with the enclosure box only!
  virtual int RegionTest(const SpaceRegion& region) const;

  // tests whether the set is void
  virtual int IsVoid() const {
    int i;
    for(i=0;i<3;i++)if(Point_min[i]>Point_max[i])return 1;
    return 0;
  }

  virtual ~SpaceRegion(){};
};



class CoordHalfspace: public SpaceRegion{
protected:
  friend class RegionSet;

  int ic, dir;
  vec_type val;

public:
  typedef CoordHalfspace* CoordHalfspaceP;
  static int cmp_func(const CoordHalfspaceP *a,const CoordHalfspaceP *b){
    if((*a)->ic==(*b)->ic){
      if((*a)->val<(*b)->val)return -1;
      else if((*a)->val>(*b)->val)return 1;
    }
    return 0;
  }

  // the sign tells which halfspace is inside:
  // 1 for point[coordn] < coordval
  // -1 otherwise
  CoordHalfspace(int coordn,vec_type coordval,int sign);

  // copy creation
  CoordHalfspace(const CoordHalfspace& other):SpaceRegion(other){
    ic=other.ic;
    dir=other.dir;
    val=other.val;
  }

  int GetDir() const {return dir;}
  vec_type GetVal() const {return val;}
  int GetCoordN() const { return ic;}

  //void SetDir(int d){dir=d;}


   // clone function
  virtual SpaceRegion *Clone() const {
    return new CoordHalfspace(*this);
  }

  virtual int Test(Vector_3* Point) const {
    if((*Point)[ic]<val){
      return dir;
    }
    else if((*Point)[ic]==val)return 0;
    else return -dir;
  }

  //virtual int RegionTest(const SpaceRegion& region);

  virtual int IsVoid() const{
    return 0;
  }
  virtual ~CoordHalfspace(){};

};



enum {SET_OR, SET_AND };

class RegionSet: public SpaceRegion{
protected:
  friend class CellRange;
  slList<CoordHalfspace *> CoordPlane[3][2];
  lList<SpaceRegion *> OtherRegs;


// functions to clean up the Set
# ifndef UNIX
# pragma argsused
# endif
  static void del_halfspc(long ind, CoordHalfspace *a){
    delete a;
  }

# ifndef UNIX
# pragma argsused
# endif
  static void del_otherreg(long ind, SpaceRegion *a){
    delete a;
  }

  virtual void test_minmax(const SpaceRegion *a, long rgbefore=-1);
  virtual void test_allminmax();

  // test function returns
  // 0 if no test is further needed
  // 1 otherwise
  // res is a result of one simple test
  // reso is a pointer to the results of preceeding tests
  // the function must accumulate the final result in *reso:
  // -2 undefined, -1 outside, 0 on the border, 2 inside
  int (*tst_func)(int *reso,int res);
  static Vector_3 *tstpoint;


  // auxilliary function to Test
  int test_list(int *reso,void *listpar) const{
    lList<SpaceRegion *> *reglist=(lList<SpaceRegion *> *)listpar;
    long curp=reglist->Rewind();
    int res, tst=1;
    //TestCount=0;
    while(!reglist->Done()){
      res=reglist->GetCur()->Test(tstpoint);
      tst=tst_func(reso,res);
      if(!tst)break;
      //TestCount++;
      reglist->Next();
    }
    reglist->SetCurPos(curp);
    return tst;
  };

  int settype;

  virtual int planeOR(int *reso) const;
  virtual int planeAND(int *reso) const;
  virtual int planeDEF(int *reso) const{
    int i,j;
    for(i=0;i<3;i++){
      for(j=0;j<2;j++){
        if(!test_list(reso,(void *)&CoordPlane[i][j]))return 0;
      }
    }
    return 1;
  }


  int nplanes;

  void init(int sttype){
    settype=sttype;
    nplanes=0;
    // default is void set
    Point_min=Vector_3(VEC_INFTY,VEC_INFTY,VEC_INFTY);
    Point_max=Vector_3(-VEC_INFTY,-VEC_INFTY,-VEC_INFTY);
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++)CoordPlane[i][j].SetFunc(&CoordHalfspace::cmp_func);
  }

  RegionSet(int sttype, int (*tstf)(int *, int)){
    init(sttype);
    tst_func=tstf;
  }

public:
  int TestCount; // to find out on which region is quit by testing

  static int tstAND(int *reso, int res){
    if(*reso==-1)return 0; // out
    if(res==-2){ // if res undefined
      *reso=-2;
      return 0;  // undefined
    }
    if(res==-1){
      *reso=-1;
      return 0; // out
    }
    /* res>=0 */
    if(*reso>0 || *reso==-2)*reso=res; // takes res
    // otherwise leaves on border
    return 1;
  }

  static int tstOR(int *reso, int res){
    if(*reso>=0)return 0; // inside or on border, no test needed
    *reso=res;
    if(res>=0)return 0;// inside or on border
    if(*reso==-2)return 0; // totally undefined
    return 1; // out
  }



  RegionSet(int sttype=SET_AND){
    init(sttype);
    switch(sttype){
      case SET_OR:
        tst_func=tstOR;
        break;
      case SET_AND:
        tst_func=tstAND;
        break;
      default:
        fatal_error("RegionSet: Undefined Set Type!\n");
    }

  }

  // copy creation
  RegionSet(const RegionSet& other):SpaceRegion(other){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++){
      CoordPlane[i][j]=other.CoordPlane[i][j];
      CoordPlane[i][j].Rewind();
      while(!CoordPlane[i][j].Done()){
        CoordPlane[i][j].SetCur(new CoordHalfspace(*CoordPlane[i][j].GetCur()));
        CoordPlane[i][j].Next();
      }
    }
    OtherRegs=other.OtherRegs;
    OtherRegs.Rewind();
    while(!OtherRegs.Done()){
      OtherRegs.SetCur(new SpaceRegion(*OtherRegs.GetCur()));
      OtherRegs.Next();
    }
    tst_func=other.tst_func;
    nplanes=other.nplanes;
    settype=other.settype;
  }

   // clone function
  virtual SpaceRegion *Clone() const {
    return new RegionSet(*this);
  }

  // adding a halfspace constraint
  // the sign tells which halfspace is inside:
  // 1 for point[coordn] < coordval
  // -1 otherwise
  // returns the Region reference pointer
  void *AddCoordPlane(int ic,vec_type coordval,int sign){
    void *ptr;
    CoordHalfspace *halfspc=new CoordHalfspace(ic,coordval,sign);
    test_minmax(halfspc);
    CoordPlane[halfspc->ic][(1+sign)/2].AppendP(halfspc,1,&ptr);
    nplanes++;
    return ptr;
  }
  void *AddCoordPlane(const CoordHalfspace *halfspc){
    void *ptr;
    test_minmax(halfspc);
    CoordPlane[halfspc->ic][(1+halfspc->dir)/2].AppendP(new CoordHalfspace(*halfspc),1,&ptr);
    nplanes++;
    return ptr;
  };

  void *AddRegion(SpaceRegion *reg, int deleg=0){
    void *ptr;
    test_minmax(reg);
    if(deleg)OtherRegs.AppendP(reg,1,&ptr);
    else OtherRegs.AppendP(reg->Clone(),1,&ptr);
    return ptr;
  }

  // removes plane by reference pointer
  // returns 1 if successfull, 0 otherwise
  int RemoveCoordPlane(void *ptr){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++){
      CoordHalfspace **plane=CoordPlane[i][j].GetByPtr(ptr);
      if(plane){
        delete *plane;
        CoordPlane[i][j].RemoveByPtr(ptr);
        nplanes--;
        test_allminmax();
        return 1;
      }
    }
    return 0;
  }

  CoordHalfspace* GetCoordPlane(void *ptr){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++){
      CoordHalfspace** plane=CoordPlane[i][j].GetByPtr(ptr);
      if(plane)return *plane;
    }
    return NULL;
  }


  int RemoveRegion(void *ptr){
    SpaceRegion **reg=OtherRegs.GetByPtr(ptr);
    if(reg){
      delete *reg;
      OtherRegs.RemoveByPtr(ptr);
      test_allminmax();
      return 1;
    }
    else return 0;
  }

  int RemoveSome(void *ptr){
    if(RemoveCoordPlane(ptr))return 1;
    else if(RemoveRegion(ptr))return 1;
    return 0;
  }

  virtual ~RegionSet(){
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++)CoordPlane[i][j].cIterate(&RegionSet::del_halfspc);
    OtherRegs.cIterate(&RegionSet::del_otherreg);
  }

  // tells whether the point is inside
  // returns 1 if inside, 0 -- on the border, -1 --outside
  virtual int Test(Vector_3* Point) const{
    if(IsVoid())return -1;
    int tstreso=-2; // undefined
    tstpoint=Point;
    switch(settype){
      case SET_OR:
        if(!planeOR(&tstreso))return tstreso;
        break;
      case SET_AND:
        if(!planeAND(&tstreso))return tstreso;
        break;
      default:
        if(!planeDEF(&tstreso))return tstreso;
        break;
    }
    test_list(&tstreso,(void *)&OtherRegs);
    return tstreso;
  }

  // deletes unnecessary planes for predefined RegionSets  (SET_OR, SET_AND)
  virtual void Refine();

  virtual long NPlanes() const{
    return nplanes;
    /*long nc=0;
    int i,j;
    for(i=0;i<3;i++)for(j=0;j<2;j++)nc+=CoordPlane[i][j].Count();
    return nc;*/
  }

  virtual long NRegs() const{
    return NPlanes()+OtherRegs.Count();
  }

  virtual SpaceRegion* GetRegion(long i){
    if(i<0 || i>=OtherRegs.Count())return NULL;
    return OtherRegs.Get(i);
  }

};



class SpaceDistribution: public RegionSet{
protected:
  SpaceRegion *reg; // limiting region for functions only
  //RegionSet *Final;

  vec_type locate_midpoint(TableFunction *distr, int n1,int n2);

  virtual SpaceDistribution* ExtractSub() const{
   SpaceDistribution *ptr=(SpaceDistribution *)Clone();
   ptr->Refine();
   return ptr;
 }


public:
  // default constructor for a homogenuous distribution inside some region
  // the region limits the definition area of the distribution
  // default is unit cube
  SpaceDistribution(const SpaceRegion* region=NULL): RegionSet(SET_AND){
    if(region){
      reg=region->Clone();
    }
    else reg=new SpaceRegion(Vector_3(0,0,0),Vector_3(1,1,1));
    AddRegion(new SpaceRegion(),1); // adding the SPACE
    //Final=new RegionSet(SET_OR);
  }

  SpaceDistribution(const SpaceDistribution& other):RegionSet(other){
    reg=other.reg->Clone();
  }

  virtual SpaceRegion* Clone() const {
    return new SpaceDistribution(*this);
  }


  SetLimitRegion(const SpaceRegion* region){
    delete reg;
    reg=region->Clone();
  }

   // getting limits for distribution functions
  Vector_3 GetFuncMin(){
    int i;
    Vector_3 a=reg->GetMin();
    for(i=0;i<3;i++)a[i]=fmax(a[i],Point_min[i]);
    return a;
  }
  Vector_3 GetFuncMax(){
    int i;
    Vector_3 a=reg->GetMax();
    for(i=0;i<3;i++)a[i]=fmin(a[i],Point_max[i]);
    return a;
  }

  virtual void *AddPlaneConstraint(int ic,vec_type val,int dir){
    return AddCoordPlane(ic,val,dir);
  }

  virtual int RemovePlaneConstraint(void *ptr){
    return RemoveCoordPlane(ptr);
  }

  virtual void *ReversePlaneConstraint(void *ptr){
    CoordHalfspace *plane=GetCoordPlane(ptr);
    if(!plane)return NULL;
    int ic=plane->GetCoordN(), dir=plane->GetDir();
    vec_type val=plane->GetVal();
    RemoveCoordPlane(ptr);
    return AddCoordPlane(ic,val,-dir);
  }


  // getting coordinate subdistributions
  # pragma argsused
  virtual void GetCoordDistr(int ic, Distribution *dens){
    // default  homogenous distribution
    TableFunction *rr=dens->get_func();
    rr->insert_begin(3);
    rr->insert_next(1.);
    rr->insert_next(1.);
    rr->insert_next(1.);
  }

  // Result is a set of resulting distributions, obtained
  // by calling ExtractSub from the constrained original one
  virtual int Subdivide(int nparts,vec_type avlength, RegionSet *Result);


  virtual ~SpaceDistribution(){
    delete reg;
  };
};





# endif

