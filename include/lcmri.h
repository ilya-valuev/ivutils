# ifndef _LCMRI_H
# define _LCMRI_H

/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2008        All Rights Reserved.
 *
 *   Author     : Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project    : ivutils
 *
 *   $Revision: 1.1 $
 *   $Date: 2008/11/11 09:10:55 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/lcmri.h,v 1.1 2008/11/11 09:10:55 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/lcmri.h,v $
$Revision: 1.1 $
$Author: valuev $
$Date: 2008/11/11 09:10:55 $
*/
/*s****************************************************************************
 * $Log: lcmri.h,v $
 * Revision 1.1  2008/11/11 09:10:55  valuev
 * added local CMRI header
 *
*******************************************************************************/

/*r @file lcmri.h 
@brief Local Cluster Multi-Range Interpolation (l-CMRI) algorithm */

class lcmri_t{
  template <class sorted_t>
  struct sort_t {
    bool operator()(const sorted_t& a, const sorted_t& b) const {
      return a.p<b.p;
    }
  };

public:
  typedef size_t iterator;
  
  struct switcher_t{
    double p;
    bool state;
    size_t ind;
    switcher_t(double p_, bool state_, size_t ind_):p(p_),state(state_),ind(ind_){}
  };
  vector<switcher_t> bonds;
  size_t mconf;

  lcmri_t():mconf(0){}

  void clear(){
    bonds.clear();
    mconf=0;
  }
  //e adds new element
  void add(double pi, bool sti){
    if(pi<0.5){
      pi=1.-pi;
      sti=!sti;
    }
    bonds.push_back(switcher_t(pi,sti,bonds.size()));
  }
  
  // gets the probability of closed bond with given initial index
  double get_closed_p(size_t i) const {
    size_t k=0;
    while(bonds[k].ind!=i)
      k++;
    return bonds[k].state ? bonds[k].p : 1.-bonds[k].p;
  }

  // gets the probability of having particular number of other existing bonds
  // presuming the given bond closed 
  // ONLY num_others=0 IS IMPLEMENTED!
  double bondconf(size_t i_closed, size_t num_others=0) const {
    if(num_others==0){
      double p=1.;
      size_t k=0, nb=nbonds();
      for(;k<nb;k++){
        if(bonds[k].ind!=i_closed){
          double pk=bonds[k].state ? bonds[k].p : 1.-bonds[k].p;
          p*=(1-pk);
        }
      }
      return p; 
    }
    else
      return 0.;
  }

  void prepare(){
    sort_t<switcher_t> pred;
    sort(bonds.begin(),bonds.end(),pred);
    size_t nmax = 8*sizeof(iterator);
    size_t nsh= min( nmax, bonds.size());
    mconf=1<<nsh;
  }

  //e returns full number of configurations (limited by sizeof(iterator))
  size_t max_conf() const {
    return mconf;
  }

  size_t nbonds() const {
    return bonds.size();
  }
  iterator begin() const{ return 0; }
  iterator end() const { return max_conf(); }

  // gets the probability and the state of the ith configuration starting from maximal
  double pconf(const iterator conf, vector<bool> &state) const {
     size_t nmax = 8*sizeof(iterator);
     nmax= min( nmax, bonds.size());
     double p=1.;
     size_t i;
     for(i=0; i< nmax; i++){
       if(conf&(1<<i)){ // flipped
         p*=1-bonds[i].p;
         state[bonds[i].ind]=!bonds[i].state;
       }
       else{ // direct
         p*=bonds[i].p;
         state[bonds[i].ind]=bonds[i].state;
       }
     }
     for(;i<bonds.size();i++){ // the rest is set as direct
       p*=bonds[i].p;
       state[bonds[i].ind]=bonds[i].state;
     }
     return p;
  }
};



# endif