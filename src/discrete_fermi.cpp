/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2014        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *   $Revision: 1.1 $
 *   $Date: 2014/07/28 08:39:49 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/discrete_fermi.cpp,v 1.1 2014/07/28 08:39:49 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/discrete_fermi.cpp,v $
$Revision: 1.1 $
$Author: valuev $
$Date: 2014/07/28 08:39:49 $
*/
/*e****************************************************************************
 * $Log: discrete_fermi.cpp,v $
 * Revision 1.1  2014/07/28 08:39:49  valuev
 * discrete fermi distributions
 *
 *
*******************************************************************************/

#include "discrete_fermi.h"

double exc_en_fermi(int nlev, int ne, double beta_de, FILE *f=NULL, double equant = .1) {
  
  double ec_i = 0.5, e_sum=0.;
  for(int i=0;i<nlev;i++){
    int level=i+1;
    /*
    int a_prev = 1.;
    for(int g=ne-1;g>=0;g--){
      double a_cur = (exp((ne-g)*beta_de) - 1. )* exp(-level*beta_de);
      a_prev = 1. - a_cur*a_prev;       
    }
    double prob = a_prev;
    
   */
    
    double prob = 0.;
    double sign = 1;
    for(int g=0;g<ne;g++){
      double prod = sign;
      
      

      for(int p=0;p<g+1;p++){
        prod *= exp((ne-p-level)*beta_de) - exp(-level*beta_de);
      }
      sign = -sign;
      prob+= prod;
    }
    e_sum+=ec_i*prob;
    if(f)
      fprintf(f,"%d %g %g %g %g\n",i,(0.5+i)*equant,prob,ec_i*prob,e_sum*equant);
    ec_i+=1.;
  }
  double e_ground = ne*ne/2.;
  return (e_sum - e_ground)*equant;
}


double exc_en_hartree(int nlev, int ne, double beta_de, FILE *f=NULL, double equant = .1) {
  double ec_i = 0.5, e_sum=0., p_sum =0.;;
  
  for(int i=0;i<nlev;i++){
    p_sum+=exp(-i*beta_de);;
  }

  for(int i=0;i<nlev;i++){  
    double prob = ne*exp(-i*beta_de)/p_sum;
    e_sum+=ec_i*prob;
    if(f)
      fprintf(f,"%d %g %g %g %g\n",i,(0.5+i)*equant,prob,ec_i*prob,e_sum*equant);
    ec_i+=1.;
  }
  double e_ground = ne/2.;
  return (e_sum - e_ground)*equant;
}


# if 0

int main(int narg, char **args){
  FILE *fl = fopen("fermi1d.dat","wt");
  fprintf(fl,"#1-num 2-Ei 3-occupation 4-Ei*occ 5-E_integ\n");
  FILE *fh = fopen("hartree1d.dat","wt");
  fprintf(fh,"#1-num 2-Ei 3-occupation 4-Ei*occ 5-E_integ\n");
  FILE *fexc = fopen("e_exc.dat","wt");
  fprintf(fexc,"#1-t 2-Eexc_fermi 3-Eexc_hartree 4-Exc_class\n");
  double equant = 11.42995;
  int ne =4, nlev =10;
  for(double t=0.2;t<2;t+=0.1){
    fprintf(fl,"#t=%g equant\n",t);
    fprintf(fh,"#t=%g equant\n",t);
    double e_exc_f =exc_en_fermi(nlev,ne,1./t,fl,equant);
    double e_exc_h =exc_en_hartree(nlev,ne,1./t,fh,equant);
    fprintf(fl,"\n");
    fprintf(fh,"\n");
    fprintf(fexc,"%g %g %g %g\n",t*equant,e_exc_f/ne,e_exc_h/ne,t*equant);
  }
  fclose(fl);
  fclose(fh);
  fclose(fexc);
  
}

#endif