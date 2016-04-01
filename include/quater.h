# include "vector_3.h"
# include "basis_3.h"

# ifndef __QUATER_H
# define __QUATER_H

class Quaternion {
 public:
  vec_type q0;
  Vector_3 vector;

  Quaternion(vec_type r0=0,vec_type r1=0,vec_type r2=0,vec_type r3=0){
   q0=r0;
   vector=Vector_3(r1,r2,r3);
  };

  Quaternion(float *vect){ q0=vect[0];vector=Vector_3(vect+1);};
  Quaternion(double *vect){ q0=vect[0];vector=Vector_3(vect+1);};
  Quaternion(vec_type a,Vector_3 vec){ q0=a;vector=vec;};
  Quaternion(Vector_3 vec){ q0=0;vector=vec;};

  vec_type& operator[](int i){
   if(!i)return vector[i-1];
   else return q0;
  };                     

  Quaternion& operator=(const Quaternion &quat){
   q0=quat.q0;
   vector=quat.vector;
   return *this;
  };                            

  Quaternion operator!() const{
    return Quaternion(q0,-vector);
  }

  Quaternion operator*(const Quaternion &q) const {
   return Quaternion(q0*q.q0-(vector*q.vector),q.vector*q0 + vector*q.q0 + vector%q.vector);
  }
  
  // rotates a vector
  Vector_3 operator()(const Vector_3 &vec) const {
    Quaternion q2=((*this)*Quaternion(vec))*(!(*this));
    return q2.vector;
  }
  /*Quaternion Zhopa(){
   Quaternion a,b,c;
   c=a*b;
   return c;
  } */
};

inline Vector_3 RotateVector(const Vector_3& vec,const Vector_3& direc,vec_type angle){
 vec_type s=sin(angle/2),c=cos(angle/2);

 Quaternion q1(c,direc*s),qv(vec),q2;
 q2=(q1*qv)*(!q1);
 return q2.vector;
}


inline Basis_3 RotateBasis(const Basis_3& B,const Vector_3& direc,vec_type angle){
 Basis_3 B2;
 int i;
 for(i=0;i<3;i++){
  B2.vec[i]=RotateVector(B.vec[i],direc,angle);
 }
 return B2;
}


void ShowRotation(Basis_3& xyz_basis,void drawfunc(void),int swappages=0);


# endif