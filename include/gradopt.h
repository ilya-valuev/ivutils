/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2012        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *****************************************************************************/
  
/*s****************************************************************************
 * $Log: gradopt.h,v $
 * Revision 1.4  2012/09/28 11:25:40  valuev
 * updated degeneracy constraints
 *
 * Revision 1.3  2012/09/20 05:38:31  valuev
 * optimizer (test version)
 *
 * Revision 1.2  2012/04/29 08:52:29  valuev
 * added power-conditional constraint
 *
 * Revision 1.1  2012/04/18 10:36:36  valuev
 * optimizer implementation (not working yet)
 *
 *
*******************************************************************************/
# ifndef GRADOPT_H
# define GRADOPT_H
  
/** @file gradopt.h 
    @brief Template classes for gradient-based optimization. */


// reformulated form optgeom.h

enum optLSEARCH_RESULT {
  optCONVERGED=1,
  optTHRESH=2,
  optTRIAL=3,
  optNODESCENT =-1,
  optINVDESCENT=-2,
  optITERLIMIT =-3,
  optPROCFAILED  =-4
};

enum optMETHODS {
  optSTEEPEST=1,
  optFLETCHER=2, //e< Fletcher-Reeves
  optPOLAK=3     //e< Polak-Ribiere
};


///\en Base optimizer class
class base_optimizer{
public:
  ///\en Set convergence creiteria:\n
  ///    max_df -- maximal function variance at convergence OR;\n
  ///    maximal position variance at convergence, \a dx is normalized by total variable number OR\n
  ///    maximal gradient norm \a max_grad at convergence.\n
  ///    Negative parameters mean that the corresponding criteria are not checked.
  ///    These conditions have to be repeated for at least \a max_convit consecutive iterations. 
  virtual int SetAccuracy(double max_df=-1, double max_dx=-1, int max_convit=2, double max_grad=-1){
    return 1;
  }

  virtual int Optimize(int maxiter=-1){
    return 1;
  }

  ///\en Function for displaying optimization status lines (legacy, avoid using this).
  ///    \param num number of symbols to display
  //virtual void show_status(int num, const char *line){
  // }
};


/// Gradient descent with parabolic minima search
template<class term_t>
class md_grad_optimizer: public base_optimizer{
public:
  typedef int (callback_t)(int res, double *x, double val, double *grad);
protected:

  int dim;
  int proc_calls;
  callback_t *pcallback;
  


  double grdval, E0, tm;
  double IncrCoeff, tStep, tStep0;
  double min_dE, min_dt;
  int method;
  int iter_lim, max_convit;
  int max_flats;
  int vout;
  FILE *f;
  
  int thresh;
  double Ethresh;
  double kdir;
  int try_shrink;
  int *fixvars;
  
protected:
  ///\en function to calcualte value
  int tFunc(double &val,double t){
    int i;
    for(i=0;i<dim;i++){
      xref[i]=x[i]-kdir*dir[i]*t;  
    }
    proc_calls++;
    term->copy_variable(0x1 /*MD_COORDS*/, false, xref.begin());
    return term->compute(&val);
  }

  ///\en makes linear search using tFunc, grdval 
  ///    @returns
  int ParabolicMin(FILE *f, int vout=1);

  typedef typename term_t::value_t value_t;
  vector<value_t> x;
  vector<value_t> xref;
  vector<value_t> grad;
  vector<value_t> grad0;
  vector<value_t> dir;
  mngptr<term_t> term;
public:


  md_grad_optimizer(term_t *term_, int managed=0):fixvars(NULL){
    init(term_,managed);
  }

  int init(term_t *term_, int managed=0){
    E0=0.;
    term.reset(term_,managed);
    if(!term.ptr())
      return 0;
    dim=(int)term->dimension();
    x.resize(dim);
    term->copy_variable(0x1 /*MD_COORDS*/,true,x.begin());
    xref=x;
    //xref.resize(dim);
    grad.resize(dim);
    grad0.resize(dim);
    dir.resize(dim);

    kdir = -1./dim; // assuming force instead of gradient

    SetVerbose(1,stdout);
    SetAccuracy(1e-5,1e-5,3);
    SetThreshold(0,0);
    SetStep(0.1,1.1);
    SetMethod(optSTEEPEST);
    max_flats=10;
    return 1;
  }


  int SetVerbose(int be_verbose, FILE *sf){
    vout=be_verbose;
    f=sf;
    return 1;
  }

  ///\en Set convergence creiteria:\n
  ///    max_dE -- maximal energy variance at convergence OR;\n
  ///    maximal position variance at convergence, \a max_dx is normalized by total variable number OR\n
  ///    maximal force norm \a max_force at convergence.\n
  ///    Negative parameters mean that the corresponding criteria are not checked.
  ///    These conditions have to be repeated for at least \a maxconvit consecutive iterations.
  virtual int SetAccuracy(double max_dE, double max_dx, int maxconvit=2, double max_force=-1.){
    min_dE=max_dE;
    min_dt=max_dx;
    max_convit=maxconvit;
    return 1;
  }

  int SetFixedVars(int *ptr){
    fixvars=ptr;
    return 1;
  }

  template <class it_t>
  int ApplyGradConstraints(it_t it){
    if(!fixvars)return 0;
    int i;
    for(i=0;i<dim;i++, ++it){
      if(fixvars[i])
        *it=0.;
    }
    return 1;
  }

  int SetThreshold(int use_thresh, double val_thresh){
    thresh=use_thresh;
    Ethresh=val_thresh;
    return 1;
  }

  int SetStep(double xstep, double coeff){
    tStep=tStep0=xstep;
    IncrCoeff=coeff;
    return 1;
  }

  int SetMethod(int meth){
    method=meth;
    return 1;
  }

  int Optimize(int maxiter=-1);
  
  //virtual int GetCurrent(double *val, double *x){
  //  *val=E0;
  //  int i;
  //  for(i=0;i<dim;i++)x[i]=xref[i];
  //  return 1;
  //}

  
};

template<class term_t>
int md_grad_optimizer<term_t>::ParabolicMin(FILE *f, int vout){
  char /*str[200],*/ str1[20];
 
  double Es,t0,t1,E1,t2,E2,Em;
  t0=0;
  t1=tStep;
  Es=E0-grdval*t1;

  
  if(vout){
    if(f){
      fprintf(f,"\nDescent:\n");
      fprintf(f,"Gradient: %g\n",grdval);
      fprintf(f,"E0=%g, Es=%g, tStep=%g\n",E0,Es,tStep);
    }
  }

  int nflat=0;
  do{
    if(tFunc(E1,t1)<0)
      return optPROCFAILED;
    // estimation of tStep
    if(fabs(Es-E1)<10*min_dE)
      tStep*=IncrCoeff;
    
    if(fabs(Es-E1)> fabs(Es-E0))
      tStep/=IncrCoeff;
    

    double coeff=1e40;
    if(fabs(E1-E0)>1e-32)
      coeff=grdval/((E1-E0)/t1);
    
    
    if(vout){
      if(f)
        fprintf(f,"Es=%g, E1=%g, tStep=%g Grd/der=%g\n",Es,E1,tStep,coeff);
    }  
    if(fabs(E1-E0)<min_dE && nflat<max_flats){
      t1=t1+tStep;;
      nflat++;
      if(vout){
        if(f)
          fprintf(f,"flat (%d of %d)\n",nflat,max_flats);
      }     
    }
    else break;
    Es=E0-grdval*t1;
  }while(1);

  int its=0;
  if(E1>E0){
    if(!try_shrink){
      if(vout){
        if(f)
          fprintf(f,"Reverse step:\n");
      }
      t2=t0-tStep;
      if(tFunc(E2,t2)<0)
        return optPROCFAILED;
      

      if(E2>E0){
        if(vout){
          if(f)
            fprintf(f,"rev E0=%g, t2=%g, E2=%g\n",E0,t2,E2);
        }
        tm=0.;
        return optNODESCENT; // recalculate command: better different direction from t0
      }
      else{
        tm=0.;
        return optINVDESCENT; // the function's descent direction is opposite to the given one
      }
    }
    else{
      if(vout){
        if(f)
          fprintf(f,"Shrinking:\n");
      }
      do{
        E2=E1;
        t2=t1;
        t1=t2/2;

        if(t1>min_dt){
          if(tFunc(E1,t1)<0)
            return optPROCFAILED;
          if(vout){
            if(f)
              fprintf(f,"sh%d E0=%g, t1=%g, E1=%g\n",its+1,E0,t1,E1);
          }
        }

        if(t1<min_dt || fabs(E1-E0)<min_dE){
          tm=t0;
          Em=E0;
          if(vout){
            if(f)
              fprintf(f,"Accuracy limit!\n");
            
            if(tFunc(E1,-tStep/4.)<0)
              return optPROCFAILED;  // Test by going to the opposite direction

            term->copy_variable(0x1 /*MD_COORDS*/, false, x.begin()); // returning to min after check
           
            if(f)
              fprintf(f,"Check point: (%g, %g)\n",(double)-tStep/4.,E1);
          }
       
          return optCONVERGED; // result is tm, Em
        }
        its++;
      }while(E1>E0);
    }
  }
  if(thresh){ // checking threshold
    if(E1<Ethresh){
      if(vout){
        if(f)
          fprintf(f,"thr! Thresh=%g, E=%g\n",Ethresh,E1);
      }
      tm=t1;
      return optTHRESH; // result (thresholded) is in t1,E1 
    }
  }
  
  t2=t1+tStep;
  if(tFunc(E2,t2)<0)
    return optPROCFAILED;
  

  double a=0.,b=0.,c=0.;
  double r0=0.,r1=0.,r2=0.;

  // now getting to descent
  its=0;
  do{
    if(its>100){ // maximal number of parabolic iterations
      if(tm>t2){
        tm=t2;
        Em=E2;
      }
      if(tm<t0){  ////?
        tm=t0;
        Em=E0;
      }
      break;
    }
    if(its){
      if(its==1)
        if(vout && f)fprintf(f,"Expanding:\n");
      t0=t1;
      E0=E1;
      t1=t2;
      E1=E2;
      t2=t1+IncrCoeff*(t1-t0);
      if(tFunc(E2,t2)<0)
        return optPROCFAILED;
    }
    if(thresh){ // checking threshold
      if(E2<Ethresh){
        //got_thresh=1;
        if(vout){
          if(f)
            fprintf(f,"thr! Thresh=%g, E=%g\n",Ethresh,E2);
        }
        tm=t2;
        return optTHRESH; // result (thresholded) is in t2, E2
      }
    }
   
    double d10=t1-t0;
    double d21=t2-t1;
    double d20=t2-t0;
    r0=1./(d10*d20);
    r1=-1./(d10*d21);
    r2=1./(d20*d21);

    a=E0*r0+E1*r1+E2*r2;
    if(a<=0){
      if(vout){
        sprintf(str1,"%d",its);
        if(f)
          fprintf(f,"%s%2s (%.g, %g) (%g,%g) (%g,%g) concave\n",(its? "ex": "VA"),(its? str1:"L "),t0,E0,t1,E1,t2,E2);
      }
      its++;
      continue;
    }
    b=-E0*r0*(t1+t2)-E1*r1*(t0+t2)-E2*r2*(t0+t1);
    tm=-b/(2*a);
    if(vout){
      sprintf(str1,"%d",its);
      if(f)
        fprintf(f,"%s%2s (%g, %g) (%g,%g) (%g,%g) %g\n",
        (its? "ex": "VA"),(its? str1:"L "),t0,E0,t1,E1,t2,E2,tm);
    }
    its++;
    
  }while(a<0 || (tm>t2 && tm-t2>2*(t2-t0)));


  c=t1*t2*r0*E0+t2*t0*r1*E1+t1*t0*r2*E2;
  Em=tm*(a*tm+b)+c;

  if(vout && f)
    fprintf(f,"tm=%g, Em=%g\n",tm,Em);

  
  //tStep=tm; // ??
  return optTRIAL; // result is in tm, Em
}


template<class term_t>
int md_grad_optimizer<term_t>::Optimize(int maxiter){
  
  char str[200];
  if(maxiter<0)
    iter_lim=1000;
  else iter_lim=maxiter;

  proc_calls = 0;


  tStep=tStep0;

  int i;
  for(i=0;i<dim;i++){ // copying coords to work with
    x[i]=xref[i];
  }
  int convit=0;
  int reset_dir=0;
  int res=0;
  double E0prev=0., gnorm=0, gnorm0=0;
  int iter=0;
  do{// through possible directions
    if(!reset_dir){

      term->copy_variable(0x1 /*MD_COORDS*/, false, xref.begin()); 
      proc_calls++;
      int res1=term->compute(&E0,dim? &grad[0] : NULL);
      if(res1<0){
        res=optPROCFAILED; // procedure call failed?
        break;
      }
      ApplyGradConstraints(grad.begin());
      gnorm=0.;
      for(i=0;i<dim;i++)
        gnorm+=grad[i]*grad[i];
      
      
      // trying to find out wheteher we converged
      if(iter>0 && ( fabs(E0-E0prev)<min_dE || fabs(tm)<min_dt)){
        convit++;
        if(convit>=max_convit)
          res=optCONVERGED; // yes, converged
      }
      else {
        convit=0;
      }
      E0prev=E0;
      if(iter>=iter_lim)
        res=optITERLIMIT;
      
      //if(pcallback)
        //pcallback(res,xref,E0,pgrad.get_ptr()); // notifying the caller about iteration step
      iter++;
      if(res==optCONVERGED|| res==optTHRESH) {
        if(vout){
          res==optCONVERGED ?
            sprintf(str,"Converged in %i iterations", proc_calls) : sprintf(str,"Threshold after %i iterations", proc_calls);
          if(f)
            fprintf(f,"%s\n",str);
        }
        break; //converged
      }
      if(iter>iter_lim)
        break; // limit reached
      if(res<0)
        break; // not converged
    }

    // descent direction
    if(iter==0 || method==optSTEEPEST || reset_dir){ // setting dir to gradient
      //pdir.copy_data(pgrad); // initial setting to direction is gradient  
      dir=grad;
      
      try_shrink=1; // this is to tell ParabolicMin to try shrinking
      if(vout && method!=optSTEEPEST && iter!=0){
        if(f)
          fprintf(f,"direction reset to gradient\n");
      }
      grdval=dim? sqrt(gnorm)/dim : 0.;
      
    }
    else{ // finding conjugate direction
      try_shrink=0; // shrinking is useless if direction is wrong
      double gamma;
      if(method==optFLETCHER){ // Fletcher-Reeves
        gamma=gnorm0/gnorm;
      }
      else{  // Polak-Ribiere
        double mul=0.;
        for(i=0;i<dim;i++)mul+=grad[i]*grad0[i];
        gamma=(gnorm0-mul)/gnorm;
      }
      grdval=0.;
      for(i=0;i<dim;i++){
        dir[i]=grad[i]+gamma*dir[i];
        grdval+=dir[i]*dir[i];
      }
      grdval=dim? sqrt(grdval)/dim : 0.;
    }
    gnorm0=gnorm;
    grad0=grad;
    res=ParabolicMin(f,vout);
    if(res>0){ // new point is accepted point
      for(i=0;i<dim;i++){ // setting the coords according to found point
        xref[i]=x[i]-kdir*tm*dir[i];
        x[i]=xref[i];
      }
    }
    else if(res<0){ // new point is rejected because it did not lead to minimization 
      for(i=0;i<dim;i++){ // copying coords back (restoring)
        xref[i]=x[i];
      }
      if(method!=optSTEEPEST && !reset_dir){ // try to reset direction
        reset_dir=1;
        continue; // give another chance
      }
    }
    reset_dir=0;
  }while(1);
  if(f && f!=stdout)
    fflush(f);
  return res;
}



# endif