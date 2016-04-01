# include <math.h>

# include "curve2d.h"
# include "common.h"

//int Curve2d::nc =0;


int Curve2d::AddPoint(float x,float y){
  input=1;
  n++;
  if(n>nlim){
    if(!reinit(n+D_NPOINTS))serror("Curve2d: AddPoint: can't increase"
				   " array size!\n");
  }
  px[n-1]=x;
  py[n-1]=y;
  if(n!=1){
    float dx=px[n-1]-px[n-2];
    float dy=py[n-1]-py[n-2];
    length+=sqrt(dx*dx+dy*dy);
  }
  return n;
}

float normalize2(float nrm,float &v1, float &v2){
  float vnorm=sqrt(v1*v1+v2*v2);
  if(fabs(vnorm)<1e-32)return 0.;
  v1/=vnorm*nrm;
  v2/=vnorm*nrm;
  return vnorm;
}

float norm2sq(float v1, float v2){
  return v1*v1+v2*v2;
}

float norm2(float v1, float v2){
  return sqrt(v1*v1+v2*v2);
}



int Curve2d::CrossPoints(float x01,float x02,float v1, float v2,
			 int **ind,float **t,float **xc, float **yc){
  if(n==0)return 0;
  Pool plind((void **)ind,(int)sizeof(int));
  Pool plt((void **)t,(int)sizeof(float));
  Pool plx((void **)xc,(int)sizeof(float));
  Pool ply((void **)yc,(int)sizeof(float));


  float fact, fact1;
  fact=(py[0]-x02)*v1-(px[0]-x01)*v2;

  
  int ncross=0;
  int i;
  for(i=1;i<n;i++){
    
    fact1=(py[i]-x02)*v1-(px[i]-x01)*v2;

    if((fact>=0 && fact1<=0) || (fact<=0 && fact1>=0)){
      // crossing point found

      float l1,l2;
      int itmp=i-1;
      plind.add(&itmp,1); //ind[ncross]=i-1;

      l1=px[i]-px[i-1];
      l2=py[i]-py[i-1];
      //normalize2(1.,l1,l2);
      
      float tmp=fact/(v2*l1-v1*l2);
     
      plt.add(&tmp,1);
      if((*t)[ncross]>1. || (*t)[ncross]<0){
	printf("Curve2d: warning: bug! t=%f tmp=%f\n",(*t)[ncross],tmp);
      }

      tmp=px[i-1]+l1*(*t)[ncross];
      plx.add(&tmp,1);
      tmp=py[i-1]+l2*(*t)[ncross];
      ply.add(&tmp,1);
      ncross++;
    }
    fact=fact1;
  }
  return ncross;
}


int Curve2d::ClosestCross(float x01,float x02,float v1, float v2,
			  int *ind,float *t,float *x, float *y){
  float *pt, *pxc, *pyc;
  int *pind;

  int ncr=CrossPoints(x01,x02,v1,v2,&pind,&pt,&pxc,&pyc);
  
  if(ncr<=0)return 0;

  int i;
  int minind=0;
  float minnrm=norm2sq(x01-pxc[0],x02-pyc[0]);
  float nrm;
  for(i=1;i<ncr;i++){
    nrm=norm2sq(x01-pxc[i],x02-pyc[i]);
    if(nrm<minnrm){
      minnrm=nrm;
      minind=i;
    }
  }
  *ind=minind;
  *t=pt[minind];
  *x=pxc[minind];
  *y=pyc[minind];
  
  free(pt);
  free(pxc);
  free(pyc);
  free(pind);

  return ncr;
}

float Curve2d::Distance(float x01,float x02,int *ind,
			float *t, float *x, float *y){
  if(n==0)return -1.; // no points
  int i;
  int minind=0;
  float minnrm=norm2sq(x01-px[0],x02-py[0]);
  float nrm;
  for(i=1;i<n;i++){
    nrm=norm2sq(x01-px[i],x02-py[i]);
    if(nrm<minnrm){
      minnrm=nrm;
      minind=i;
    }
  }

  /*
  *t=0;
  *x=px[minind];
  *y=py[minind];
  *ind=minind;
  return sqrt(minnrm);*/
  
  float da=-1, db=-1;
  float va1=0,va2=0,vb1=0,vb2=0,prja=0,prjb=0,nva,nvb;
  float r1,r2,nr;


  r1=x01-px[minind];
  r2=x02-py[minind];
  nr=norm2(r1,r2);
  if(nr<1e-32){
    *ind=minind;
    *t=0;
    *x=px[minind];
    *y=py[minind];
    return nr;
  }


  if(minind>0){
    va1=px[minind-1]-px[minind];
    va2=py[minind-1]-py[minind];
    nva=norm2(va1,va2);
    
    if(nva<1e-32)da=nr;
    else{
      prja=(va1*r1+va2*r2)/nr/nva;
      if(prja>0){
	va1=prja*va1*nr/nva-r1;
	va2=prja*va2*nr/nva-r2;
	da=norm2(va1,va2);
      }
      else da=-1;
    }
  }
  if(minind<n-1){
    vb1=px[minind+1]-px[minind];
    vb2=py[minind+1]-py[minind];
    nvb=norm2(vb1,vb2);
    
    if(nvb<1e-32)db=nr;
    else{
      prjb=(vb1*r1+vb2*r2)/nr/nvb;
      if(prjb>0){
	vb1=prjb*vb1*nr/nvb-r1;
	vb2=prjb*vb2*nr/nvb-r2;
	db=norm2(vb1,vb2);
      }
      else db=-1;
    }
  }

    

  if(da<0 && db<0){ // all angles too big
    *ind=minind;
    *t=0;
    *x=px[minind];
    *y=py[minind];
    return nr;
  }

  if(da<0 || (db>0 && db<=da)){
    *ind=minind;
    *t=prjb;
    *x=px[minind]+vb1;
    *y=py[minind]+vb2;
    return db;
  }

  *ind=minind-1;
  *t=1-prja;
  *x=px[minind]+va1;
  *y=py[minind]+va2;

  return da; 
}

int BoxedCurve::AddPoint(float x, float y){
  if(Band<0){
    crv_out=0;
    added=1;
    //printf("hoho\n");
    return Curve2d::AddPoint(x,y);
  }
  if(!IN_FRAME(x,Xmin-Band,Xmax+Band) || !IN_FRAME(y,Ymin-Band,Ymax+Band)){
    if(!crv_out){
     crv_out=1;
     added=1;
     return Curve2d::AddPoint(x,y);
    }
    else{
      crv_out=2;
      xlast=x;
      ylast=y;
      added=0;
      return n;
    }
  }
  if(crv_out==2)AddPoint(xlast,ylast);
  crv_out=0;
  added=1;
  return Curve2d::AddPoint(x,y);
}













