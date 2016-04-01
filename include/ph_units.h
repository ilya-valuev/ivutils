# ifndef PH_UNITS_H
# define PH_UNITS_H


# include "common.h"
  


class dim_t{
  static char names[];  //="MLTQ";
  static char sstr[250];
public:
  int x[8];

  void check_d(int u, int d){
     if(d<=0)serror("Dimension: incorrect dimension set: %d/%d !\n",u,d);
  }

public:
  dim_t(){
    x[0]=0;
    x[1]=1;
    x[2]=0;
    x[3]=1;
    x[4]=0;
    x[5]=1;
    x[6]=0;
    x[7]=1;
  }
  void M(int u,int d){
    check_d(u,d);
    x[0]=u;
    x[1]=d;
  }
  void L(int u,int d){
    check_d(u,d);
    x[2]=u;
    x[3]=d;
  }
  void T(int u,int d){
    check_d(u,d);
    x[4]=u;
    x[5]=d;
  }
  void Q(int u,int d){
    check_d(u,d);
    x[6]=u;
    x[7]=d;
  }

  int operator!=(dim_t &dim);
  int operator==(dim_t &dim){
    return !((*this)!=dim);
  }
  char *str();
  dim_t &operator=(dim_t d);
  friend dim_t operator-(dim_t d1, dim_t d2);
  friend dim_t operator+(dim_t d1, dim_t d2);
  dim_t operator-();
};


class ph_unit{

protected:
  double value;
  dim_t dim;
  int strict;

  virtual void set_dim(){
    dim.M(0,1);
    dim.L(0,1);
    dim.T(0,1);
    dim.Q(0,1);
    strict=0;
    printf("ncall\n");
  }

 

public:
  
  ph_unit(double val=0){
    value=val;
    set_dim();
    //printf("Initialized:\n%d/%d\n%d/%d\n\n",
    //	   dim.x[0],dim.x[1],dim.x[2],dim.x[3]);
  }

  ph_unit(ph_unit& U){
    value=U.value;
    dim=U.dim;
    strict=U.strict;
  }

  ph_unit &operator=(ph_unit U){
    if(strict){
      if(dim!=U.dim){
	eprintf("ph_unit: trying to initialize unit with dimension: %s\n",
		dim.str());
	serror("by unit with dimension: %s\n",U.dim.str());
      }
    }
    value=U.value;
    dim=U.dim;
    // strict=U.strict;  ??
  }

  friend ph_unit operator+(ph_unit , ph_unit );
  friend ph_unit operator-(ph_unit , ph_unit );
  friend ph_unit operator*(double , ph_unit);
  friend ph_unit operator*(ph_unit, double);
  friend ph_unit operator/(double , ph_unit );
  friend ph_unit operator/(ph_unit, double);
  friend ph_unit operator*(ph_unit U1, ph_unit U2);
  friend ph_unit operator/(ph_unit U1, ph_unit U2);
};


class ph_length:public ph_unit{
  virtual void set_dim(){
    dim.M(0,1);
    dim.L(1,1);
    dim.T(0,1);
    dim.Q(0,1);
    //printf("vcall\n");
    strict=1;
  }
public:
  ph_length(double val=0.){
    value=val;
    set_dim();
  }
  ph_length(ph_unit U){
    set_dim();
    printf("eql %d\n",this->strict);
    ph_unit *tmp=this;
    *tmp=U;
  }
  /*
  ph_unit &operator=(ph_unit U){
    ph_unit &tmp=*this;
    tmp=U;
    return tmp;
  } */
};



# endif 





