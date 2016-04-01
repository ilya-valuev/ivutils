# include <string.h>
# include <stdio.h>

# include "ph_units.h"


int dim_t::operator!=(dim_t &U){
  int i;
  for(i=0;i<8;i++)if(x[i]!=U.x[i]){
    //printf("ne: %s %s\n",str(),U.str());
    return 1;
  }
  return 0;
}

char *dim_t::str(){
  int i;
  char tmp[30];
  sstr[0]=0;
  for(i=0;i<4;i++){
    if(x[2*i]!=0){
      if(x[2*i+1]!=1)sprintf(tmp,"[%c]^(%d/%d)",names[i],
				x[2*i],x[2*i+1]);
      else sprintf(tmp,"[%c]^%d",names[i],x[2*i]);
      strcat(sstr,tmp);
    }
  }
  if(sstr[0]==0)sprintf(sstr,"0");
  return sstr;
}

dim_t & dim_t::operator=(dim_t d){
  int i;
  for(i=0;i<8;i++)x[i]=d.x[i];
  return *this;
}

dim_t operator+(dim_t d1, dim_t d2){
  dim_t dr;
  int i,div;
  for(i=0;i<4;i++){
    dr.x[2*i]=d1.x[2*i]*d2.x[2*i+1]+d1.x[2*i+1]*d2.x[2*i];
    dr.x[2*i+1]=d1.x[2*i+1]*d2.x[2*i+1];
    div=com_div(dr.x[2*i+1],dr.x[2*i]);
    dr.x[2*i]/=div;
    dr.x[2*i+1]/=div;
  }
  return dr;
}

dim_t operator-(dim_t d1, dim_t d2){
  dim_t dr;
  int i,div;
  for(i=0;i<4;i++){
    dr.x[2*i]=d1.x[2*i]*d2.x[2*i+1]-d1.x[2*i+1]*d2.x[2*i];
    dr.x[2*i+1]=d1.x[2*i+1]*d2.x[2*i+1];
    div=com_div(dr.x[2*i+1],dr.x[2*i]);
    dr.x[2*i]/=div;
    dr.x[2*i+1]/=div;
  }
  return dr;
}

dim_t dim_t::operator-(){
  int i;
  for(i=0;i<4;i++){
    x[2*i]=-x[2*i];
  }
  return *this;
}

void dim_t::mul(int u, int d){
  int div=com_div(u,d);
  u/=div;
  d/=div;
  int i;
  for(i=0;i<4;i++){
    x[2*i]*=u;
    x[2*i+1]*=d;
  }
}


char dim_t::names[]="MLTQ";
char dim_t::sstr[250];


void check_lin_dim(dim_t &d1,dim_t &d2){
  if(d1!=d2){
    eprintf("ph_unit: trying to perform linear operation with\n"
	    "units of different dimension:\n%s and ",d1.str());
    serror("%s\n",d2.str());
  }
}


ph_unit operator+(ph_unit U1, ph_unit U2){
  check_lin_dim(U1.dim,U2.dim);
  ph_unit Ur;
  Ur.value=U1.value+U2.value;
  Ur.dim=U1.dim;
  return Ur;
}


ph_unit operator-(ph_unit U1, ph_unit U2){
  check_lin_dim(U1.dim,U2.dim);
  ph_unit Ur;
  Ur.value=U1.value-U2.value;
  Ur.dim=U1.dim;
  return Ur;
}

ph_unit operator*(double a, ph_unit U){
  U.value*=a;
  return U;
}

ph_unit operator*(ph_unit U, double a){
  return a*U;
}

ph_unit operator/(ph_unit U,double a){
  U.value=U.value/a;
  return U;
}

ph_unit operator/(double a, ph_unit U){
  U.value=a/U.value;
  U.dim=-U.dim;
  return U;
}

ph_unit operator*(ph_unit U1, ph_unit U2){
  ph_unit Ur;
  Ur.value=U1.value*U2.value;
  Ur.dim=U1.dim+U2.dim;
  return Ur;
}


ph_unit operator/(ph_unit U1, ph_unit U2){
  ph_unit Ur;
  Ur.value=U1.value/U2.value;
  Ur.dim=U1.dim-U2.dim;
  return Ur;
}
  

ph_unit ph_pow(ph_unit U,int u,int d=1){
  ph_unit Ur;
  double r=((double)u)/d;
  Ur.value=pow(U.value,r);
  Ur.dim=U.dim.mul(u,d);
  return Ur;
}



const ph_length ph_m(1.), ph_cm(1e-2), ph_mm(1e-3), ph_A(1e-10);


void main(){
  printf("main\n");
  ph_unit A, B;
  //ph_length L;
  //L=2.3*ph_m;
  ph_length L=2.3*ph_m,L1=2*ph_cm;
  printf("step1:\n");

  //L=L*L;
  
  L1=L/2+L;
  A=B;
  L1=ph_cm;
  B=L;
  L=L1;
  L=3*L+L1*5;
  printf("Now:\n");
  L=(L*L)/(L1*L*L);
}








