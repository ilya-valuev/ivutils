//# include <iostream.h>
# include <stdio.h>
# include <math.h>

# include "vector_n.h"

Vector_Np& Vector_Np::operator=(float *vect){
  int i;
  for(i=0;i<size;i++)v[i]=vect[i];
  return *this;
}

Vector_Np& Vector_Np::operator=(double *vect){
  int i;
  for(i=0;i<size;i++)v[i]=vect[i];
  return *this;
}

Vector_Np& Vector_Np::operator=(Vector_Np& vect){
  chksz(vect.size);
  int i;
  for(i=0;i<size;i++)v[i]=vect.v[i];
  return *this;
};


Vector_Np Vector_Np::operator+(Vector_Np& vect){
  chksz(vect.size);
  int i;
  Vector_Np result(size);
  for(i=0;i<size;i++)result.v[i]=v[i]+vect.v[i];
  return result;
};

Vector_Np Vector_Np::operator-(Vector_Np& vect){
  chksz(vect.size);
  int i;
  Vector_Np result(size);
  for(i=0;i<size;i++){
   result.v[i]=v[i]-vect.v[i];
  }
  return result;
};

vec_typen Vector_Np::operator*(Vector_Np& vect){
  chksz(vect.size);
  int i;
  vec_typen result=0.;
  for(i=0;i<size;i++)result+=v[i]*vect.v[i];
  return result;
};

Vector_Np Vector_Np::operator*(vec_typen coeff){
  int i;
  Vector_Np result(size);
  for(i=0;i<size;i++)result.v[i]=v[i]*coeff;
  return result;
};

Vector_Np Vector_Np::operator/(vec_typen coeff){
  int i;
  Vector_Np result(size);
  for(i=0;i<size;i++)result.v[i]=v[i]/coeff;
  return result;
};

Vector_Np operator*(vec_typen coeff, Vector_Np& vect){
  return vect*coeff;
};

Vector_Np Vector_Np::operator-(){
  Vector_Np r(size);
  int i;
  for(i=0;i<size;i++)r.v[i]=-v[i];
  return r;
}

Vector_Np& Vector_Np::operator*=(vec_typen coeff){
  int i;
  for(i=0;i<size;i++)v[i]*=coeff;
  return *this;
}

Vector_Np& Vector_Np::operator/=(vec_typen coeff){
  int i;
  for(i=0;i<size;i++)v[i]/=coeff;
  return *this;
}

Vector_Np& Vector_Np::operator+=(Vector_Np& vec){
  int i;
  for(i=0;i<size;i++)v[i]+=vec.v[i];
  return *this;
}

Vector_Np& Vector_Np::operator-=(Vector_Np& vec){
  int i;
  for(i=0;i<size;i++)v[i]-=vec.v[i];
  return *this;
}



vec_typen Vector_Np::norm(){
 vec_typen sum=0.;
 int i;
 for(i=0;i<size;i++)sum+=v[i]*v[i];
 return sqrt(sum);
}

vec_typen Vector_Np::normalize(){
 vec_typen norm;
 norm=this->norm();
 if(norm!=0)(*this)/=norm;
 return norm;
}

                                                                             

void Vector_Np::print(char *form){
 char *acform="%f ";
 if(form)acform=form;
 int i;
 for(i=0;i<size;i++){
   printf(acform,v[i]);
 }
 printf("\n");
}


// multiplies symmetric band matrix of with n1-1 subdiagonals by vec
int BandMul(Vector_Np &resvec,int ns,vec_typen *band_mat,Vector_Np &vec){
  int i,j;

  resvec.chksz(vec.size);

  for(j=0;j<vec.size;j++){
    resvec.v[j]=0;
    int i1=fmax(0,j-ns+1);
    int i2=fmin(j+ns-1,vec.size-1);

    for(i=i1;i<j;i++){
      resvec.v[j]+=band_mat[i*ns+j-i]*vec.v[i];
    }
    for(i=j;i<=i2;i++){
      resvec.v[j]+=band_mat[j*ns+i-j]*vec.v[i];
    }
  }
  return 1;
}







