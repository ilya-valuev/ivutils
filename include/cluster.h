# ifndef CLUSTER_H
# define CLUSTER_H



class pair_dist{
public:   
  pair_dist *next;
  double r;
  int c1;
  int c2;
  
  pair_dist(int i1,int i2,double d){
    c1=i1;
    c2=i2;
    r=d;
    next=NULL;
  }
  //void* operator new(size_t sz,int i1,int i2,double d);
  void *operator new(size_t sz);
  void operator delete(void *buf);
};


typedef pair_dist *pair_distP;


class ClusterFunc{
  int search_place(pair_dist* &op,int i1, int i2);
  void insert_distance(int c1,int c2,double r);
  int update_distance(int c1,int c2,double r);
  int join_distance(int resting);
  pair_dist *entry;
public:
  int *ind;
  int n;
  int n_clust;
  double distance;

  ClusterFunc(int num,double dist);
  ClusterFunc(const ClusterFunc &other);
  ~ClusterFunc();

  // returns the distance between clusters
  // i1, i2 -- corresponding point indicies
  // -1. -- if not known
  // -2. -- if at least one cluster does not exist
  // -3. -- internal bug
  double r(int cl1,int cl2, int &i1, int &i2);

  //e returns the maximum distance between clusters
  //e 0 -- if there is 1 cluster or less
  double rmax(int &i1,int &i2);

  //e informs about the new interpoint distance
  //e @return the number of clusters
  int inform(int i, int j, double d, double comp_dist=-1);
  int nchain();
  void reset();
};
# endif




