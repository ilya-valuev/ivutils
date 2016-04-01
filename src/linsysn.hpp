# include<math.h>
# include<stdlib.h>
# include<complex>
# include"linsysn.h"     
//# include"common.h"

using namespace std;


template<class T>
int linsys_cpp(T **matr,T *vect,int num,T *x){
//int linsys_cpp(T matr[][3],T *vect,int num,T *x){
 int *ind=new int[num];
 if(!ind)return -1; /* mememory allocation failure */
 for(int i=0;i<num;i++)ind[i]=i;

 for(int i=0;i<num;i++){
  x[i]=0;
  int nm=i;
  double max=abs(matr[i][ind[i]]);
  for(int k=i+1;k<num;k++){
    double f=abs(matr[i][ind[k]]);
    if(f>max){
       max=f;
       nm=k;
    }
  }
  /*if(max<1e-20 && )return 0;*/
  /*else{*/
  if(nm!=i){
    int ti=ind[i];
    ind[i]=ind[nm];
    ind[nm]=ti;
    T t=vect[i];
    vect[i]=vect[nm];
    vect[nm]=t;
//   swap_var(ind[i],ind[nm],int);
//   swap_var(vect[i],vect[nm],T);
  }
  if(max>1e-20){
   T mul=1./matr[i][ind[i]];
   for(int k=i+1;k<num;k++){
    if(abs(matr[i][ind[k]])!=0){
     vect[k]-=vect[i]*matr[i][ind[k]]*mul;
     for(int j=i+1;j<num;j++)matr[j][ind[k]]-=matr[j][ind[i]]*matr[i][ind[k]]*mul;
    }
    matr[i][ind[k]]=0;
   }
  }
  else matr[i][ind[i]]=0;
 }
 for(int i=num-1;i>=0;i--){
  if(abs(matr[i][ind[i]])<1e-20){
   if(abs(vect[i])<1e-20)x[i]=vect[i]=1;
   else return 0;
  }
  else x[i]=vect[i]=vect[i]/matr[i][ind[i]];
  matr[i][ind[i]]=1;
  for(int j=0;j<i;j++){
   if(abs(matr[i][ind[j]])!=0){
    vect[j]-=matr[i][ind[j]]*vect[i];
    matr[i][ind[j]]=0;
   }
  }
 }
 delete[]ind;
 return 1;
}

