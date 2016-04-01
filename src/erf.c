# include <math.h>
# include <stdio.h>

# define SQRTPI  1.772454
# define ZERO    0
# define ONE     1
# define TWO     2
# define THREE   3
# define TWENTY  20
# define HALF    0.5

# define SIGN(x,y) ( (y)>=0 ? (x) : -(x) )

double erf(double x);

double erf(double x){
 int ncfc=8, ncfd=8,j;
 double XUP=4.0, xv, x2, bjp2, bjp1, bj,er;
 double C[8]={                  \
  1.944907,4.2019E-2,-1.8687E-2,5.129E-3,-1.068E-3,\
  1.74E-4,-2.1E-5,2.0E-6  \
 };
 double D[8]={              \
  1.483110,-3.01071E-1,6.8995E-2,-1.3916E-2,2.421E-3,\
 -3.66E-4,4.9E-5,-6.0E-6 \
 };


 xv = fabs(x);
 if(xv < XUP){

  if( xv > TWO){
    x2 = TWO - TWENTY/(xv+THREE);

/*     SUMMATION  */
    bjp2 = ZERO;
    bjp1 = C[ncfc-1];

    for(j = ncfc - 1;;j--){
     bj = x2*bjp1 - bjp2 + C[j-1];
     if(j==1)break;
     bjp2 = bjp1;
     bjp1 = bj;
    }

    x2 =( HALF*(bj - bjp2)/xv)*exp(-x*x)/SQRTPI;
    er = (ONE-x2)*SIGN(ONE,x);
  }
  else{
   x2 = x*x - TWO;
/*     SUMMATION   */
   bjp2 = ZERO;
   bjp1 = D[ncfd-1];
   for(j = ncfd - 1;;j--){
    bj = x2*bjp1 - bjp2 + D[j-1];
    if(j==1)break;
    bjp2 = bjp1;
    bjp1 = bj;
   }
   er = HALF*(bj-bjp2)*x;
  }
 }
 else er = SIGN(ONE,x);
 return er;
}


double erfc(double x){
 return 1.-erf(x);
}

