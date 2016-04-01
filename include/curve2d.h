# ifndef CURVE2D_H
# define CURVE2D_H

# include "common.h"

# define D_NPOINTS   20 

// class requires initialization !

class Curve2d{
protected:
  TableFunction fx, fy;
  int input;
  virtual int end_input(){
    if(!input)return 0;
    input=0;
    fx=TableFunction(n,px);
    fy=TableFunction(n,py);
    return 1;
  }
public:
  int nlim;
  int n;
  float *px;
  float *py;
  float length;
  // static int nc;

  Curve2d(){
   
    //printf("con0: %d ",nc++);
    px=py=NULL;
    n=nlim=0;
    length=0.;
    input=1;
  }

  virtual ~Curve2d(){
    //printf("des: %d\n",nc--);
    if(px)free(px);
    if(py)free(py);
  }
  
  int init(int num){
    //printf("ini ");
    if(px)free(px);
    if(py)free(py);

    input=1;
    length=0.;
    nlim=num;
    n=0;
    px= (float *)malloc(num*sizeof(float));
    py= (float *)malloc(num*sizeof(float));
    if(!px || ! py)return 0;
    else return 1;
  }


  int reinit(int num){
    if(nlim==0)return init(num);
    px= (float *)realloc(px, num*sizeof(float));
    py= (float *)realloc(py, num*sizeof(float));

    if(!px || !py)return 0;
    nlim=num;
    return 1;
  }

  Curve2d(int num){
    if(!init(num))serror("Curve2d: memory allocation error.\n");
  } 

  virtual int AddPoint(float x,float y);

  //int GetMidPoint(int ip,float t,float *x,float *y){
  //  if(t<0 || l> length)return -1;

  // returns the number of cross points and their specs
  // with the line (x01,x02)+(v1,v2)*tau
  int CrossPoints(float x01,float x02,float v1, float v2,
		  int **ind,float **t,float **x, float **y);

  // returns the number of cross points and
  // the closest one in ind, t ,x, y
  int ClosestCross(float x01,float x02,float v1, float v2,
		   int *ind,float *t,float *x, float *y);

  // returns the closest distance between the curve and 
  // point (x01,y01)
  // in ind, t, x, y -- curve reference of "perpendicular"
  float Distance(float x01, float x02,int *ind,float *t, float *x, float *y);

  float Dx(int ind, float t){
    if(ind<0 || ind>=n)return 0;
    end_input();
    return fx.Dy(((float)ind+t)/n);
  }
  
  float Dy(int ind, float t){
    if(ind<0 || ind>=n)return 0;
    end_input();
    return fy.Dy(((float)ind+t)/n);
  }

  float Dl(int ind, float t){
    float dx=Dx(ind,t);
    float dy=Dy(ind,t);
    return sqrt(dx*dx+dy*dy);
  }

  float X(int ind, float t){
    if(ind<0 || ind>=n)return 0;
    end_input();
    return fx((float)ind+t);
  } 

  float Y(int ind, float t){
    if(ind<0 || ind>=n)return 0;
    end_input();
    return fy((float)ind+t);
  } 
 
};



class BoxedCurve: public Curve2d{
  int crv_out;
  float xlast, ylast;
public:
  float Xmin, Xmax, Ymin, Ymax;
  float Band;

  int added;

  BoxedCurve():Curve2d(),Xmin(0),Xmax(0),Ymin(0),Ymax(0),Band(-1){}

  int init(int num,float x1,float y1, float x2, float y2, float bw){

    Xmin=x1;
    Ymin=y1;
    Xmax=x2;
    Ymax=y2;
    Band=bw;
    if(Xmin>Xmax){
      x1=Xmin;
      Xmin=Xmax;
      Xmax=x1;
    }
    if(Ymin>Ymax){
      y1=Ymin;
      Ymin=Ymax;
      Ymax=y1;
    }
    crv_out=1;
    added=0;
    return Curve2d::init(num);
  }

  BoxedCurve(int num,float x1,float y1, float x2, float y2, float bw){
    init(num,x1,y1,x2,y2,bw);
  }

  BoxedCurve(int num){
    crv_out=0;
    added=0;
    Band=-1;
    Curve2d::init(num);
  }
   

  virtual int AddPoint(float x,float y);

  //virtual ~BoxedCurve(){}
};


float normalize2(float nrm,float &v1, float &v2);


# endif
















