# define _USE_MATH_DEFINES
# include <cmath>
# include "splines.h"
# include "utiltl.h"


double f1(size_t n, double x){
  double sgn= n%4>1 ? -1 : 1;
  double f= n%2 ? cos(x) : sin(x);
  return f*sgn;
}

/// f2=f1(r), r=sqrt(x*x+y*y)
double f2(size_t nx, double x, size_t ny, double y){
  double r=sqrt(x*x+y*y);
  if(nx==0 && ny==0)
    return f1(0,r);
  if(nx+ny==1)
    return f1(1,r)*(nx== 1 ? x: y)/r;
  else if(nx==1 && ny==1)
    return (f1(2,r)-f1(1,r)/r)*x*y/(r*r);
  else if(nx+ny==2)
    return (f1(2,r)*(nx==2 ? x*x : y*y) +f1(1,r)*(nx==2 ? y*y : x*x)/r)/(r*r);
  return 0.;
}

# define TEST_1D  0
# define TEST_1D_AKIMA 1
# define TEST_1DM 0
# define TEST_2D  0
# define TEST_2DM 0


int main(){
  FILE *f;
  GridSpline<double> spl(1,10);

# if TEST_1D
  //-------------- random 1D grid -------------------------
  // starting build 
  spl.BeginBuild();

  // setting random x values
  for(int i=0;i<10;i++){
    rand();
    rand();
    spl.set_arg(0,i,2*M_PI*(1.-(double)(rand())/RAND_MAX));
  }
  spl.CommitGrdAxis(0); // commit arguments
  
  f=fopen("spl1or.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl''(x) 9-f''(x)\n");
  for(int i=0;i<10;i++){
    double xx=spl.get_arg(0,i);
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,f1(0,xx),f1(0,xx),f1(1,xx),f1(1,xx),f1(2,xx),f1(2,xx),f1(3,xx),f1(3,xx));
    spl.grd_value(i)=f1(0,xx);
  }
  fclose(f);
  spl.EndBuild();
  // default ends (cubic)
  f=fopen("spl1r.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl''(x) 9-f''(x)\n");
  for(double xx =-4., dx=0.01;xx<2*M_PI+4.;xx+=dx){
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,spl.der(0,xx),f1(0,xx),spl.der(1,xx),f1(1,xx),spl.der(2,xx),f1(2,xx),spl.der(3,xx),f1(3,xx));
    
  }
  fclose(f);
  
  //------------- regular 1D grid ------------------------
  spl.BeginBuild();
  // quadratic ends
  spl.SetExtrapolation(0,0,1);
  spl.SetAxis(0,0.,2*M_PI);
  spl.SetEnds(0,1,3);
  f=fopen("spl1o.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl'''(x) 9-f'''(x)\n");
  for(int i=0;i<10;i++){
    double xx=2*M_PI*i/9;
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,f1(0,xx),f1(0,xx),f1(1,xx),f1(1,xx),f1(2,xx),f1(2,xx),f1(3,xx),f1(3,xx));
    spl.grd_value(i)=f1(0,xx);
  }
  fclose(f);
  spl.EndBuild();
  f=fopen("spl1.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl''(x) 9-f''(x)\n");
  for(double xx =-4., dx=0.01;xx<2*M_PI+4.;xx+=dx){
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,spl.der(0,xx),f1(0,xx),spl.der(1,xx),f1(1,xx),spl.der(2,xx),f1(2,xx),spl.der(3,xx),f1(3,xx));
  }
  fclose(f);

# endif // TEST_1D
  
# if TEST_1D_AKIMA
  //-------------- irregular 1D grid -------------------------
  double table[]={8,20, 15,20, 23,20, 31,10, 40,40, 46,10, 55,0, 63,0, 64,0, 65,0 };

  // starting build 
  spl.BeginBuild();
  spl.SetAxis(0,table[0],table[9*2]);
  spl.SetEnds(0,1,3);
  spl.SetInterpolation(0,1);

  // setting x values
  for(int i=0;i<10;i++){
    spl.set_arg(0,i,table[2*i]);
    spl.grd_value(i)=table[2*i+1];
  }
  spl.CommitGrdAxis(0); // commit arguments
  
  f=fopen("spl1ak_data.d","wt");
  fprintf(f,"#1-x 2-data\n");
  for(int i=0;i<10;i++){
    fprintf(f,"%g %g \n",table[2*i],table[2*i+1]);
  }
  fclose(f);
  spl.SetMode(SPLINE_NATURAL|SPLINE_AKIMA);
  spl.EndBuild();
  // default ends (cubic)
  
  f=fopen("spl1ak.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-spl'(x) 4-spl''(x) 5-spl'''(x)\n");
  for(double xx =0., dx=0.1;xx<66.;xx+=dx){
    fprintf(f,"%g %g %g %g %g\n",xx,spl.der(0,xx),spl.der(1,xx),spl.der(2,xx),spl.der(3,xx));
  }
  fclose(f);
  
# endif // TEST_1D_AKIMA






  int np=30;
  vector<double> px, py;

# if TEST_1DM
  spl.BeginBuild(0);
  // quadratic ends
  spl.SetExtrapolation(0,0,1);
  spl.SetAxis(0,0.,2*M_PI);
  spl.SetEnds(0,1,3);
  for(int i=0;i<np;i++){
    px.push_back(2*M_PI*(1.-(double)(rand())/RAND_MAX));
  }
  sort(px.begin(),px.end());

  int ind=0;
  spl.BeginSubgrid(&ind,1);
  f=fopen("spl1om.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl'''(x) 9-f'''(x)\n");
  for(int i=0;i<np;i++){
    double xx=px[i]; //2*M_PI*(i+0.3)/9;
    double cc=1;
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,f1(0,xx),f1(0,xx),f1(1,xx),f1(1,xx),f1(2,xx),f1(2,xx),f1(3,xx),f1(3,xx));
    //spl.grd_value(i)=f1(0,xx);
    spl.AddEquation(1,&xx,&cc,f1(0,xx));
  }
  fclose(f);
  spl.EndSubgrid();
  spl.EndBuild();
  f=fopen("spl1m.d","wt");
  fprintf(f,"#1-x 2-spl(x) 3-f(x) 4-spl'(x) 5-f'(x) 6-spl''(x) 7-f''(x) 8-spl''(x) 9-f''(x)\n");
  for(double xx =-4., dx=0.01;xx<2*M_PI+4.;xx+=dx){
    fprintf(f,"%g %g %g %g %g %g %g %g %g\n",xx,spl.der(0,xx),f1(0,xx),spl.der(1,xx),f1(1,xx),spl.der(2,xx),f1(2,xx),spl.der(3,xx),f1(3,xx));
  }
  fclose(f);

# endif // TEST_1DM


  //----- regular 2D grid -----------------------
  Tensor<double> grd2(2,10,10);
  double extra=2.,dd=0.2; 
# if TEST_2D
  // setting 10x10 grid values
  f=fopen("spl2o.d","wt");
  fprintf(f,"#1-x 2-y 3-f00 4-f10 5-f01 6-f11 7-f20 8-f02\n");
  for(int i=0;i<10;i++){
    double xx=2*M_PI*i/9;
    for(int j=0;j<10;j++){
      double yy=2*M_PI*j/9;
      double f00=f2(0,xx,0,yy);
      double f10=f2(1,xx,0,yy);
      double f01=f2(0,xx,1,yy);
      double f11=f2(1,xx,1,yy);
      double f20=f2(2,xx,0,yy);
      double f02=f2(0,xx,2,yy);
      fprintf(f,"%g %g %g %g %g %g %g %g\n",xx,yy,f00,f10,f01,f11,f20,f02);
      grd2(i,j)=f00;
    }
    fprintf(f,"\n");
  }
  fclose(f);

  // linking grid to spline
  spl.SetFuncGrid(&grd2);
  // starting the build: the flag is to specify axis build order (x first)
  spl.BeginBuild(0x1);
  // x and y  axes
  spl.SetAxis(0,0.,2*M_PI);
  spl.SetAxis(1,0.,2*M_PI);
  // pbc along x and y
  spl.SetEnds(0,2,2);
  spl.SetEnds(1,2,2);
  spl.EndBuild();

  // testing spline
  spl.SetExtrapolation(0,3,3);
  spl.SetExtrapolation(1,3,3);
 
  f=fopen("spl2r.d","wt");
  fprintf(f,"#1-x 2-y 3-s00 4-s10 5-s01 6-s11 7-s20 8-s02 9-f00 10-f10 11-f01 12-f11 13-f20 14-f02\n");
  for(double xx=0-extra; xx<2*M_PI+extra;xx+=dd){
    for(double yy=0-extra; yy<2*M_PI+extra;yy+=dd){
      double s00=spl.der(0,xx,0,yy);
      double s10=spl.der(1,xx,0,yy);
      double s01=spl.der(0,xx,1,yy);
      double s11=spl.der(1,xx,1,yy);
      double s20=spl.der(2,xx,0,yy);
      double s02=spl.der(0,xx,2,yy);

      double f00=f2(0,xx,0,yy);
      double f10=f2(1,xx,0,yy);
      double f01=f2(0,xx,1,yy);
      double f11=f2(1,xx,1,yy);
      double f20=f2(2,xx,0,yy);
      double f02=f2(0,xx,2,yy);
      fprintf(f,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",xx,yy,s00,s10,s01,s11,s20,s02,f00,f10,f01,f11,f20,f02);
    }
    fprintf(f,"\n");
  }
  fclose(f);

# endif // TEST_2D

# if TEST_2DM
  // linking grid to spline
  spl.SetFuncGrid(&grd2);
  // starting the build
  spl.BeginBuild(0);
  // x and y  axes
  spl.SetAxis(0,0.,2*M_PI);
  spl.SetAxis(1,0.,2*M_PI);
  // pbc along x and y
  spl.SetEnds(0,2,2);
  spl.SetEnds(1,2,2);

  
  //----- random 2D grid -----------------------
  np=200;
  for(int i=0;i<np;i++){
    px.push_back(2*M_PI*(1.-(double)(rand())/RAND_MAX));
    py.push_back(2*M_PI*(1.-(double)(rand())/RAND_MAX));
  }

  /*np=100;
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++){
      px.push_back(2*M_PI*i/9);
      py.push_back(2*M_PI*j/9);
    }
  }*/

  int ind=0;
  spl.BeginSubgrid(&ind,1);
  f=fopen("spl2om.d","wt");
  fprintf(f,"#1-x 2-y 3-f00 4-f10 5-f01 6-f11 7-f20 8-f02\n");
  for(int i=0;i<np;i++){
    double xx=px[i];
    double yy=py[i]; //2*M_PI*j/9;
    double f00=f2(0,xx,0,yy);
    double f10=f2(1,xx,0,yy);
    double f01=f2(0,xx,1,yy);
    double f11=f2(1,xx,1,yy);
    double f20=f2(2,xx,0,yy);
    double f02=f2(0,xx,2,yy);
    fprintf(f,"%g %g %g %g %g %g %g %g\n",xx,yy,f00,f10,f01,f11,f20,f02);
    
    double xa[]={xx,yy}, cc=1;
    spl.AddEquation(1,xa,&cc,f00);
    //grd2(i,j)=f00;
    
    //fprintf(f,"\n");
  }
  fclose(f);
  spl.EndSubgrid();
  spl.EndBuild();

  // testing spline
  spl.SetExtrapolation(0,3,3);
  spl.SetExtrapolation(1,3,3);
  
  f=fopen("spl2rm.d","wt");
  fprintf(f,"#1-x 2-y 3-s00 4-s10 5-s01 6-s11 7-s20 8-s02 9-f00 10-f10 11-f01 12-f11 13-f20 14-f02\n");
  for(double xx=0-extra; xx<2*M_PI+extra;xx+=dd){
    for(double yy=0-extra; yy<2*M_PI+extra;yy+=dd){
      double s00=spl.der(0,xx,0,yy);
      double s10=spl.der(1,xx,0,yy);
      double s01=spl.der(0,xx,1,yy);
      double s11=spl.der(1,xx,1,yy);
      double s20=spl.der(2,xx,0,yy);
      double s02=spl.der(0,xx,2,yy);

      double f00=f2(0,xx,0,yy);
      double f10=f2(1,xx,0,yy);
      double f01=f2(0,xx,1,yy);
      double f11=f2(1,xx,1,yy);
      double f20=f2(2,xx,0,yy);
      double f02=f2(0,xx,2,yy);
      fprintf(f,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",xx,yy,s00,s10,s01,s11,s20,s02,f00,f10,f01,f11,f20,f02);
    }
    fprintf(f,"\n");
  }
  fclose(f);

# endif // TEST_2DM


  return 0;
}