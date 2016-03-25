/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.9 $
 *   $Date: 2012/09/07 14:26:06 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/seqpack.cpp,v 1.9 2012/09/07 14:26:06 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/seqpack.cpp,v $
$Revision: 1.9 $
$Author: valuev $
$Date: 2012/09/07 14:26:06 $
*/
/*s****************************************************************************
 * $Log: seqpack.cpp,v $
 * Revision 1.9  2012/09/07 14:26:06  valuev
 * removed using namespace from h files
 *
 * Revision 1.8  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.7  2010/06/12 17:57:31  valuev
 * some workflow coding
 *
 * Revision 1.11  2009/02/27 00:40:58  lesha
 * 2 calles put_record_end bug is fixed
 *
 * Revision 1.10  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.5  2008/04/09 17:31:25  valuev
 * Added bands and charges, modified wequil analyser
 *
 * Revision 1.4  2008/02/27 13:49:47  valuev
 * corrected variable delink from namespace
 *
 * Revision 1.3  2008/02/21 14:02:51  valuev
 * Added parametric methods
 *
 * Revision 1.2  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.7  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.6  2007/02/20 10:11:32  valuev
 * Debugging with g++ 2.96
 *
 * Revision 1.5  2007/02/19 23:50:09  valuev
 * New makefile style; got rid of some compiler warnings
 *
 * Revision 1.4  2006/12/03 22:46:22  lesha
 * Error in put_record_end is corrected
 *
 * Revision 1.3  2006/12/02 01:53:43  lesha
 * Error in put_record_end is corrected
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/

# include <string>
# include "stdio.h"
# include "seqpack.h"

using namespace std;

/// extracts the next value at a given position 
/// pos must be less than or equal ind.size()-2
/// @returns the value or -1 if end of sequence or the postion is incorrect
/// modifies pos (the next reading position) and incrcount
int get_next_record(const vector<int> &ind, int &pos, int &incrcount, int single){
  do{
    int ic=pos;
    int cval0=ind[pos];
    //int *pind=(int *)(&ind[0]);
    if(cval0>=0){
      ic++;
      int cval=ind[ic];
      if(cval<0){
        if(incrcount<-cval){ 
          incrcount++;
        }
        else{  
          if(ind[ic+1]<0)pos=ic+2;
          else pos=ic+1;
          incrcount=0;
          continue;
        }
        if(ind[ic+1]<0){
          cval0+=(incrcount-1)*(-ind[ic+1]-1);
        }
        else{
          cval0+=(incrcount-1)*single;
        }
      }
      else{
        pos=ic;
        incrcount=0;
      }
      return cval0;
    }
    else{
      return -1; // end of record or incorrect position
    }
  }while(1);
  return -1;
}

/// specifies the next value to be recorded
/// only nonnegative values are possible
/// @retuns    0 if the value was not packed (just stored)
///            1 packed with selected increment
///            2 packed with other increment
///            3 packed as continued sequence
///           -1 tried to store after end-of sequence 
///           -2 tried to store negative value
int put_next_record(vector<int> &ind, int cur, int single){
  if(cur<0)return -2; //negative, not allowed
  int last=(int)ind.size()-1;
  if(last<=0){ // this is the first or the second element
    ind.push_back(cur);
  }
  else{
    if(ind[last]<0){ // incremented
      if(ind[last-1]<0){// the increment is not <single>
        if(ind[last-2]<0)return -1; // end of sequence ?
        int inc0=-ind[last]-1;
        int lval=ind[last-2]+inc0*(-ind[last-1]-1);
        if(cur-lval!=inc0)ind.push_back(cur); // not confroming to the increment
        else{
          --(ind[last-1]); // one more to the count
          return 3;
        }
      }
      else{ // the increment is <single>
        int lval=ind[last-1]+(-ind[last]-1)*single;
        if(cur-lval!=single)ind.push_back(cur); // not confroming to the increment
        else{
          --(ind[last]); // one more to the count
          return 3;
        }
      }
    }
    else{ // not incremented
      if(ind[last-1]<0){ // only one value in the buffer
        ind.push_back(cur); // wait till the next value comes
      }
      else{
        int inc0=ind[last]-ind[last-1];
        int inc=cur-ind[last];
        if(inc0==inc){ // increments coinside 
          if(inc==single){ // the increment is <single>, packing count only
            ind[last]=-3;
            return 1;
          }
          else if(inc>=0){ // only if they are non negative 
            // the increment is not <single>, packing increment and count
            ind[last]=-3;
            ind.push_back(-(inc+1));
            return 2;
          }
        }
        ind.push_back(cur); // store unpacked
      }
    }
  }
  return 1;
}

/// puts the end of record indicator 
int put_record_end(vector<int> &ind,int single){
  int last=(int)ind.size()-1;
  if (last>=0) {
    if (ind[last]<0){ // alrady in increment sequence
      if(last>0){
        if(ind[last-1]>=0){ // this was a single increment
          // modyfing it to double-increment
          ind.push_back(-(single+1));
        }
        else if(last>1 && ind[last-2]<0){
          return 1; // put_record_end called earlier
        }
      }
      else {
        if(last==0)
          return 1; // put_record_end called earlier
        ind.push_back(-1); // !!! this is impossible !
      }
      ind.push_back(-1);
    }
    else{
      ind.push_back(-1);
      ind.push_back(-1);
      ind.push_back(-1);
    }
  }
  else {
    ind.push_back(-1); // if the vector is empty
  }
  return 1;
}


/// function to test integer sequence packing
int test_seqpack(int len,  double seqprob, int seqlenav, int spec_val, double specprob){
  vector<int> unpacked(len);
  int_pack packed(spec_val);

  // generating the random sequence
  int i=0;
  // starting value
  int vcur=12, interval, seq=0;
  packed.start_record();
  do{
    if(i!=0){
      if(seq<=0){
        // sequence or not
        seq=(((double)rand())/RAND_MAX < seqprob ? 1 : 0);
        if(seq){
          // the length of the sequence
          seq=(int)(seqlenav*(double)rand()/RAND_MAX);
        }
        // special interval or not
        int spec=(((double)rand())/RAND_MAX < specprob ? 1 : 0);
        if(!spec){
          interval=(int)(100.*((double)rand())/RAND_MAX);
        }
        else interval=spec_val;
      }
      else seq--;
      int v=vcur+interval;
      if(v<0)break; // end of integer range, can not be packed further
      vcur=v;
    }
    unpacked[i]=vcur;
    if(i==16329){
      i++;
      i--;
    }
    packed.next_record(vcur);
    i++;
  }while(i<len);
  packed.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",i,(int)packed.size(),int(packed.packed_size()/sizeof(int)));
  printf("Comparing:\n");
  int_pack::iterator it=packed.begin(), end=packed.end();
  int res=1;
  for(i=0;i<(int)unpacked.size();i++){
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    int vu=unpacked[i];
    int vp=*it;
    //printf("%d: %d  %d\n",i,vu,vp);
    if(vu!=vp){
      printf("! %d: %d  %d\n",i,vu,vp);
      res=0;
    }
    if(i==16328){
      i++;
      i--;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  return res;
}


int add_ranges(int_pack &packer, const string &spec, int imin, int imax, bool finalize,const string &delims){
  int res=1;
  if(spec=="all"){
    for(int i=imin;i<=imax;i++)
      packer.next_record(i);
  }
  else{
    if(delims.size()<2)
      return -1;
    char comma[]=",", range[]="-";
    comma[0]=delims[0];
    range[0]=delims[1];
    char *buff=new char[spec.size()+1];
    strncpy(buff,spec.c_str(),spec.size()+1);
    
    
    char *str1=strtok(buff,comma);
    while(str1){
      int n1, n2;
      if(sscanf(str1,"%d",&n1)!=1){
        res= -2;
        break;
      }
      char *p=strstr(str1,range);
      if(p){
        *p=0;
        if(sscanf(p+1,"%d",&n2)!=1){
          res=-2;
          break;
        }
        *p=range[0];
      }
      else
        n2=n1;
      if(n1<0 || n2<0){
        res= -3;
        break;
      }
      int dn= n1==n2? 0 : n1<n2 ? 1 : -1; // this is to be replaced by reading optional step
      n2+=dn;
      do{
        if(imin<=n1 && n1<=imax) // range control
          packer.next_record(n1);
        n1+=dn;
      }while(n1!=n2);
      str1=strtok(NULL,comma);
    }
    delete [] buff;
  }
  if(res>0 && finalize)
    packer.end_record();
  return res; 
}
