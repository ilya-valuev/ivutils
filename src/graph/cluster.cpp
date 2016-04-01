# include <stdlib.h>
# include "common.h"
# include "cluster.h"


/*
void* pair_dist::operator new(int i1,int i2,double r){
  pair_dist *tmp;
  tmp=(pair_dist *)malloc(sizeof(pair_dist));
  if(!tmp)serror("pair_dist: Error allocating pair.\n");
  tmp->next=NULL;
  tmp->r=r;
  tmp->c1=i1;
  tmp->c2=i2;
  return (void *)tmp;
}*/


int num=0;

void* pair_dist::operator new(size_t sz){
  pair_dist *tmp;
  tmp=(pair_dist *)malloc(sz);
  if(!tmp)serror("pair_dist: Error allocating pair.\n");
  tmp->next=NULL;
  num++;
  //printf("N:%d ",num);
  return (void *)tmp;
}


void pair_dist::operator delete(void *buff){
  if(!buff)serror("pair_dist::delete : Deleting zero pointer!\n");
  free(buff);
  num--;
  //printf("N:%d ",num);
}



ClusterFunc::ClusterFunc(int num,double dist){
  distance=dist;
  n=num;
  n_clust=0;
  ind= new int[n];
  int i;
  for(i=0;i<n;i++)ind[i]=-1;
  entry=new pair_dist(1,-1,0.);
  entry->next=NULL;
}

ClusterFunc::ClusterFunc(const ClusterFunc &other){
  distance=other.distance;
  n=other.n;
  n_clust=other.n_clust;
  ind= new int[n];
  int i;
  for(i=0;i<n;i++)ind[i]=other.ind[i];
  pair_dist *oe=other.entry, *pe=NULL;
  entry=NULL;
  while(oe){
    entry=new pair_dist(oe->c1,oe->c2,oe->r);
    if(pe)pe->next=entry;
    pe=entry;
    oe=oe->next;
  }
  if(entry)entry->next=NULL;
}

void ClusterFunc::reset(){
  pair_dist *prev=entry->next, *cur;

  while(prev){
    cur=prev;
    prev=prev->next;
    delete cur;
  }
  int i;
  for(i=0;i<n;i++)ind[i]=-1;
  entry->next=NULL;
  n_clust=0;
}




ClusterFunc::~ClusterFunc(){
  reset();
  delete entry;
  delete [] ind;
}

int ClusterFunc::nchain(){
  pair_dist *prev=entry;
  int n=0;

  while(prev){
    prev=prev->next;
    n++;
  }
  return n;
}

// returns the distance between clusters
// i1, i2 -- corresponding point indicies
// -1. -- if not known
// -2. -- if at least one cluster does not exist
// -3. -- internal bug
double ClusterFunc::r(int cl1,int cl2, int &i1, int &i2){
  if(cl1<0 || cl1>=n_clust)return -2.;
  if(cl2<0 || cl2>=n_clust)return -2.;

  pair_dist *prev=entry;
  if(!prev)return -3; // bug

  pair_dist *cur=prev->next;
  while(cur){
    if((ind[cur->c1]==cl1 && ind[cur->c2]==cl2 ) ||
       (ind[cur->c2]==cl1 && ind[cur->c1]==cl2 )){
      i1=cur->c1;
      i2=cur->c2;
      return cur->r;
    }
    cur=cur->next;
  }
  return -1.;
}


// returns the maximum distance between clusters
// 0 -- if there is 1 cluster or less
double ClusterFunc::rmax(int &i1,int &i2){
  double rm=0; 
  pair_dist *cur=entry->next;

  while(cur){
    //printf("(%d, %d) %f\n",cur->c1,cur->c2,cur->r);
    if(rm< cur->r){
      rm=cur->r;
      i1=cur->c1;
      i2=cur->c2;
    }
    cur=cur->next;
  }
  return rm;
}

  
int ClusterFunc::search_place(pair_dist* &op,int i1, int i2){
  op=entry;
  pair_dist *p=entry->next;
  while(p!=NULL){
    if(p->c1>i1)return 0;
    if(p->c1==i1){
      if(p->c2==i2)return 1;
      if(p->c2>i2)return 0;
    }
    op=p;
    p=op->next;
  }
  return 0;
}
      


void ClusterFunc::insert_distance(int c1,int c2,double r){
  //printf("i (%d %d) %f\n",c1,c2,r);
  int i1=fmin(c1,c2);
  int i2=fmax(c1,c2);
  pair_dist *c_dist= new pair_dist(i1,i2,r);
  
  pair_dist *pl;
  search_place(pl,i1,i2);
  c_dist->next=pl->next;
  pl->next=c_dist; 
}

int ClusterFunc::update_distance(int c1,int c2,double r){
  //printf("u (%d %d) %f\n",c1,c2,r);
  if(ind[c1]==ind[c2])return 0; // same cluster
  pair_dist *prev=entry;
  if(!prev)return 0; // bug

  pair_dist *cur=prev->next;
  while(cur){
    if((ind[cur->c1]==ind[c1] && ind[cur->c2]==ind[c2] ) ||
       (ind[cur->c2]==ind[c1] && ind[cur->c1]==ind[c2] )){
      if(cur->r>r){ // substituting
	//printf("c1:%d ",nchain());
	prev->next=cur->next; // deleting existent
	delete cur;
        //printf("c2:%d ",nchain());
	insert_distance(c1,c2,r); // adding new
        //printf("c3:%d ",nchain());
	return 1;
      }
      return 0; // not needed
    }
    prev=cur;
    cur=cur->next;
  }
  // adding to pool
  insert_distance(c1,c2,r);
  
  return -1;
}

int ClusterFunc::join_distance(int resting){
  //printf("j %d\n",resting);
  pair_dist *prev=entry, *cur=prev->next;
  
  double *rmin= new double[n_clust];
  pair_dist **rmentr= new pair_distP[n_clust];
  if(!rmentr || !rmin)serror("Update_distance: MAE\n");
  int i;
  for(i=0;i<n_clust;i++)rmin[i]=-1.;

  int dt=0;
  
  while(cur){
    //printf("je\n");
    if(cur->r>0){ // skip being deleted
      if(ind[cur->c1]==resting){
	//printf("j1\n");
	if(ind[cur->c2]==resting){
	  //printf("j2\n");
	  //prev->next=cur->next; //two belong to the joined cluster
	  //delete cur;
	  //cur=prev;
	  //printf("del ind %d (%d, %d) %f\n",resting,cur->c1,cur->c2,cur->r);
	  cur->r=-1;
	  dt++;
	}
	else{  // only first belong to the joined cluster
	  //printf("j3\n");
	  if(rmin[ind[cur->c2]]<0){
	    //printf("j4\n");
	    rmin[ind[cur->c2]]=cur->r;
	    rmentr[ind[cur->c2]]=prev;
	  }
	  else if(rmin[ind[cur->c2]]>cur->r){// replacing min
	   /* printf("del min %d-%d (%d, %d) %f\n",resting,ind[cur->c2],
		 rmentr[ind[cur->c2]]->next->c1,
		 rmentr[ind[cur->c2]]->next->c2,
		 rmentr[ind[cur->c2]]->next->r);*/

	    rmentr[ind[cur->c2]]->next->r=-1;
	    dt++;
	    rmentr[ind[cur->c2]]=prev;
	    //printf("j5\n");
	    // tmp=rmentr[ind[cur->c2]]->next;
	    //printf("j5.1\n"); 
	    //rmentr[ind[cur->c2]]->next=tmp->next;
	    //printf("j5.2\n");
	    rmin[ind[cur->c2]]=cur->r;
	    //printf("j5.3\n");
	    //rmentr[ind[cur->c2]]=prev;
	    //printf("j5.4\n");
	    //delete tmp;
	    //printf("j5.5\n");

	  }
	  else{  // deleting cur
	    //printf("j6\n");
	    //printf("del min %d-%d (%d, %d) %f\n",resting,ind[cur->c2],cur->c1,cur->c2,cur->r);
	    cur->r=-1;
	    dt++;
	    //prev->next=cur->next;
	    //delete cur;
	    //cur=prev;

	  }
	}
      }
      else if(ind[cur->c2]==resting){ // only second belong to the joined cluster
       //printf("j7\n");
       if(rmin[ind[cur->c1]]<0){
         //printf("j8\n");
         rmin[ind[cur->c1]]=cur->r;
         rmentr[ind[cur->c1]]=prev;
       }
       else if(rmin[ind[cur->c1]]>cur->r){// replacing min
         // printf("j9\n");
         /*printf("del min %d-%d (%d, %d) %f\n",resting,ind[cur->c1],
            rmentr[ind[cur->c1]]->next->c1,
            rmentr[ind[cur->c1]]->next->c2,
            rmentr[ind[cur->c1]]->next->r);*/
         rmentr[ind[cur->c1]]->next->r=-1;
             dt++;
         rmin[ind[cur->c1]]=cur->r;
         //tmp=rmentr[ind[cur->c1]]->next;
         //rmentr[ind[cur->c1]]->next=tmp->next;
         //rmin[ind[cur->c1]]=cur->r;
         rmentr[ind[cur->c1]]=prev;

         //delete tmp;
       }
       else{  // deleting cur
         //printf("del min %d-%d (%d, %d) %f\n",resting,ind[cur->c1],cur->c1,cur->c2,cur->r);
         cur->r=-1;
         dt++;
         //printf("j10\n");
         //prev->next=cur->next;
         //delete cur;
         //cur=prev;
       }
      }
    }
    prev=cur;
    cur=prev->next;
    //printf("j11\n");
  }

  int del=0;
  prev=entry;
  cur=prev->next;
  while(cur){
    if(cur->r<0){ // deleting
      prev->next=cur->next;
      delete cur;
      del++;
      cur=prev;
    }
    prev=cur;
    cur=cur->next;
  }
  // printf("Bal:%d-%d ",dt,del);

  delete[] rmin;
  delete[] rmentr;
  //printf("j12\n");
  return 1;
}


int ClusterFunc::inform(int i, int j, double d, double comp_dist){
  if(comp_dist<0)comp_dist=distance;
  //printf("bc:%d ",nchain());
  int same;
  if(d<comp_dist)same=1;
  else same=0;

  if(ind[i]<0){
    if(ind[j]<0){
      ind[i]=n_clust; // making new cluster
      if(!same)n_clust++;// making the second if not together
      ind[j]=n_clust++;
    }
    else{
      if(same)ind[i]=ind[j]; // attaching to old
      else{
	ind[i]=n_clust++; // making new cluster
      }
    }
    if(!same)insert_distance(i,j,d);
  }
  else{
    if(ind[j]<0){
      if(same)ind[j]=ind[i]; // attaching to old
      else{
	ind[j]=n_clust++; // making new cluster
	insert_distance(i,j,d);
      }
    }
    else{
      if(same && ind[i]!=ind[j]){ // joining two clusters
        int deleted=fmax(ind[i],ind[j]);
        int resting=fmin(ind[i],ind[j]);
        int k;
        for(k=0;k<n;k++){
          if(ind[k]<0)continue;
          if(ind[k]==deleted)ind[k]=resting;
          else{
            if(ind[k]>deleted)ind[k]--;
          }
        }
        n_clust--;
        join_distance(resting);
      }
      if(!same)update_distance(i,j,d);
    }
  }
  return n_clust;
}











