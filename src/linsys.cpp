# include<math.h>
# include<stdlib.h>
# include"linsys.h"     
# include"common.h"
     
/* vectors are in columns of matr */
/* returns -1 if allocation error, 0 if no solution, 1 if OK*/

int linsys(lstype **matr,lstype *vect,int num,lstype *x){   
 int i,j,k,nm;
 lstype mul,f,max;
 int *ind=(int *)malloc(num*sizeof(int));
 if(!ind)return -1; /* mememory allocation failure */
 for(i=0;i<num;i++)ind[i]=i;


 for(i=0;i<num;i++){
  x[i]=0.;
  nm=i;
  max=fabs(matr[i][ind[i]]);
  for(k=i+1;k<num;k++)if( (f=fabs(matr[i][ind[k]])) > max){
   max=f;
   nm=k;
  }
  /*if(max<1e-20 && )return 0;*/
  /*else{*/
  if(nm!=i){
   swap_var(ind[i],ind[nm],int);
   swap_var(vect[i],vect[nm],lstype);
  }
  if(max>1e-20){
   mul=1./matr[i][ind[i]];
   for(k=i+1;k<num;k++){
    if(matr[i][ind[k]]!=0.){
     vect[k]-=vect[i]*matr[i][ind[k]]*mul;
     for(j=i+1;j<num;j++)matr[j][ind[k]]-=matr[j][ind[i]]*matr[i][ind[k]]*mul;
    }
    matr[i][ind[k]]=0.;
   }
  }
  else matr[i][ind[i]]=0.;
 }
 for(i=num-1;i>=0;i--){
  if(fabs(matr[i][ind[i]])<1e-20){
   if(fabs(vect[i])<1e-20)x[i]=vect[i]=1.;
   else return 0;
  }
  else x[i]=vect[i]=vect[i]/matr[i][ind[i]];
  matr[i][ind[i]]=1.;
  for(j=0;j<i;j++){
   if(matr[i][ind[j]]!=0.){
    vect[j]-=matr[i][ind[j]]*vect[i];
    matr[i][ind[j]]=0.;
   }
  }
 }
 free(ind);
 return 1;
}
