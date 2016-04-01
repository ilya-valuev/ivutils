/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.9 $
 *   $Date: 2012/06/29 10:50:12 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/listiter.h,v 1.9 2012/06/29 10:50:12 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/listiter.h,v $
$Revision: 1.9 $
$Author: valuev $
$Date: 2012/06/29 10:50:12 $
*/
/*e****************************************************************************
 * $Log: listiter.h,v $
 * Revision 1.9  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.4  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.7  2006/12/20 14:29:33  valuev
 * Updated workflow, sync with FDTD
 *
 * Revision 1.3  2006/10/27 20:41:01  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.6  2006/10/14 05:10:15  valuev
 * adapted for vs9 compilation
 *
 * Revision 1.5  2006/10/09 12:58:56  bogomolov
 * patch for compilling under unix
 *
 * Revision 1.4  2006/07/21 16:22:03  valuev
 * Added Tight Binding for graphite+O
 *
 * Revision 1.3  2006/07/14 15:50:52  valuev
 * Added functional mdPotentialSet
 *
 * Revision 1.2  2006/03/14 10:32:17  valuev
 * Added SetControls support for many components,
 * improved tcpenfine, added GRASP interface
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:11  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/
# ifndef __LISTITER_H
# define __LISTITER_H

# include <stdlib.h>
# include "common.h"

// template classes for iterable List

template <class Item>
class basicList {
protected:
  long curpos;
public:
  basicList(){
    curpos=0;
  };

  basicList(const basicList& other){
    curpos=other.curpos;
  }

  virtual basicList *Clone() const =0;

  virtual long Count() const=0;
  virtual Item& Get(long index) const=0;

  // converting lists
  virtual void SetList(const basicList* other){
    basicList *temp;
    int rm=0;
    if(other==this){
      temp=other->Clone();
      rm=1;
    }
    else temp=(basicList *)other;
    Remove_all();
    long curp=temp->Rewind();
    while(!temp->Done()){
      Append(temp->GetCur());
      temp->Next();
    }
    if(rm)delete temp;
    else temp->SetCurPos(curp);
  }


  // returns index
  // appends na items to the end of the list
  virtual long Append(const Item &, long na=1)=0;

  // removes na items starting from index
  // if na<0 removes up to the end
  // returns new count
  virtual long Remove(long index, long na=1)=0;

  virtual void Remove_all()=0;


  // iterations
  virtual long GetCurPos() const{ return curpos;}

  // returns curpos at the call moment
  virtual long Rewind(){
    long curp=curpos;
    curpos=0;
    return curp;
  }

  virtual inline int Done() const{
    if(curpos<Count())return 0;
    return 1;
  }

  // returns the current curpos and shifts 1 step forwards if possible
  virtual long Next(){
    long curp=curpos;
    if(curp<Count())curpos++;
    return curp;
  }

  virtual void SetCurPos(long curp){
    if(curp<0 || curp>Count())fatal_error("basicList: current list position %ld "
                                         "out of range [0-%ld]!\n",curp,Count());
    curpos=curp;
  }

  virtual Item GetCur() const =0;
  virtual void SetCur(Item )=0;

  // const Iterate, does not change the list
  // stops iterations if func returns 0
  virtual void cbIterate(int (*func)(long,Item)){
    long curp=Rewind();
    while(!Done()){
      if(!func(GetCurPos(),GetCur()))break;
      Next();
    }
    SetCurPos(curp);
  }

  virtual int IsThere(const Item& si){
    int ret=0;
    long curp=Rewind();
    while(!Done()){
      Item a=GetCur();
      if(a==si){
        ret=1;
        break;
      }
      Next();
    }
    SetCurPos(curp);
    return ret;
  }


  // const Iterate, does not change the list
  virtual void cIterate(void (*func)(long,Item)){
    long curp=Rewind();
    while(!Done()){
      func(GetCurPos(),GetCur());
      Next();
    }
    SetCurPos(curp);
  }

  // non-const iterate
  virtual void Iterate(Item (*func)(long,Item)){
    //if(!Count())fatal_error("basicList: trying to iterate void list!\n");
    long curp=Rewind();
    while(!Done()){
      SetCur(func(GetCurPos(),GetCur()));
      Next();
    }
    SetCurPos(curp);
  }
};

/*
template <class Item>
del_contents(basicList<Item>){

};*/


//e one-way-linked list of Items
//e if Items are pointers
//e the Items are not copied, pointers must be valid
template <class Item>
class lList: public basicList<Item> {
  using basicList<Item>::curpos;
public:
  class lPair {
  public:
    Item it;
    lPair *next;
  public:
    lPair *Next(){
      return next;
    }
    Item& GetData(){
      return it;
    }
  } *links;
  long count;
protected:
  //e returns the item previous to the one indexed by index
  //e the search is started from the element poined by first having index start_i
  //e start_i must be  <= index
  //e returns NULL if index is out of range [1,count-1]
  lPair *locate_prev(long index, lPair *first=NULL, long start_i=0) const;

  // last element
  lPair *last;

  void init(){
    count=0;
    last=links=NULL;
    curelm=links;
  }

  // iterations variables
  lPair *curelm;
  void rewind_it(){
    if( basicList<Item>::curpos<0)basicList<Item>::curpos=0;
    if(curpos>count)curpos=count;

    if(curpos==0)curelm=links;
    else if(curpos!=count)curelm=locate_prev(curpos)->next;
    else if(curpos-1>0)curelm=locate_prev(curpos-1)->next;
    else curelm=links;
  }

  void equate(const lList& other);

  // for profiling purposes
  //int rr;
  
public:
  /*lPair* create_lPair(){
    rr++;
    return new lPair;
  }*/

  // constructors
  lList(){
    //rr=0;
    init();
  }


  lList(const lList& other):basicList<Item>(other){
    equate(other);
  }
  lList& operator=(const lList &other){
    if(this!=&other){
      Remove_all();
      equate(other);
    }
    return *this;
  }

  virtual basicList<Item>* Clone() const{
    return new lList(*this);
  }


  virtual ~lList(){
    Remove_all();
  }

  virtual void Remove_all(){
    lPair *elm=links, *elmn;
    while(elm){
      elmn=elm->next;
      delete elm;
      elm=elmn;
    }
    count=curpos=0;
    last=links=curelm=NULL;
  }

  // general from basicList
  virtual inline long Count() const { return count; }

  // for speed increse
  inline int Done() const{
    if(curpos<count)return 0;
    return 1;
  }

  virtual Item& Get(long index) const;

  // appends na items to the end of the list
  virtual long Append(const Item &a, long na=1){
    return AppendP(a,na);
  }

  // returns also a pointer to the newly appended item
  virtual long AppendP(const Item &a, long na=1, void **appelm=NULL);

  virtual int RemoveByPtr(void *ptr){
    lPair *elm=links, *prev=NULL;
    long index=0;
    while(elm){
      if(ptr==(void *)elm){ // removing
        if(prev)prev->next=elm->next;
        else links=elm->next;

        if(ptr==last){
          if(prev)last=prev;
          else last=links;
        }
        delete elm;
        count--;
        if(curpos>=index)rewind_it();
        return 1;
      }
      prev=elm;
      elm=elm->next;
      index++;
    }
    return 0;
  }

  // gets item by pointer and returns index if needed
  virtual Item* GetByPtr(void *ptr, long *indq=NULL) const{
    lPair *elm=links, *prev=NULL;
    long index=0;
    while(elm){
      if(ptr==(void *)elm){ // found
        if(indq)*indq=index;
        return &elm->it;
      }
      prev=elm;
      elm=elm->next;
      index++;
    }
    return NULL;
  }
  virtual long Remove(long index, long na=1);

  // this needs to be optimized
  /*virtual int RemoveCur(){
    if(!curelm)return 0;
    if(!prevelm){
      links=curelm->next;
    }

  }*/

  // iterations

  // returns curpos at the call moment
  virtual long Rewind(){
    curelm=links;
    return basicList<Item>::Rewind();
  }

  //e returns the current curpos and shifts 1 step forwards if possible
  virtual inline long Next(){
    long curp=curpos;
    if(curelm)if(curelm->next){
      curelm=curelm->next;
    }  
    if(curpos<count)curpos++;
    return curp;
  }

  //e removes the entry at curpos so that the next one is at curpos
  //e if the end is reached curpos will be decreasing, deleting end element each time
  //e returns the current curpos 
  virtual inline long RemoveCur(){
    if(count){
      Remove(curpos);
    }
    return curpos;
  }

  virtual void SetCurPos(long curp){
    if(curpos==curp)return;

    long oldp=curpos;
    basicList<Item>::SetCurPos(curp);
    if(oldp<curpos){
      if(curpos<count){
        if(curpos==0)curelm=links;
        else curelm=locate_prev(curpos,curelm,oldp)->next;
      }
      else{
        if(curpos-1==0)curelm=links;
        else if(curpos==count){
          curelm=last;
        }
        else{
          curelm=locate_prev(curpos-1,curelm,oldp);
          if(!curelm)fatal_error("lList.SetCurPos: Invalid position %ld, list count %ld!\n",
            curp,count);
          curelm=curelm->next;
        }
      }
    }
    else rewind_it();
  }

  virtual Item GetCur() const{
    return curelm->it;
  };

  virtual Item GetLast() const{
    if(last)return last->it;
    else return Item(0);
  }

  virtual lPair *GetCurNode() const {
    return curelm;
  };


  virtual void SetCur(Item a){
    curelm->it=a;
  };

};



// sorted list
//typedef template <class Item> int (*cmpf)(const Item *,const Item *);


template <class Item>
class slList: public lList<Item> {
  typedef int (*cmpf)(const Item *,const Item *);
  //int (*cmp_func)(const Item *,const Item *);
  cmpf cmp_func;

//------------
using lList<Item>::Remove_all;
using lList<Item>::SetCurPos;
using lList<Item>::GetCurPos;
using lList<Item>::Remove;
public:
  typedef typename lList<Item>::lPair lPair;
  /*class lPair {
    public:
    Item it;
    lPair *next;
    public:
    lPair *Next(){
      return next;
    }
    Item& GetData(){
      return it;
    }
  };*/
private:  
//------------

# ifdef _BCC
# pragma argsused
# endif
  static int default_cmp(const Item *a,const Item *b){
    return 0;
  }

  /*
  lPair *topelm;
  virtual void locate_topelm(){
    if(count<=1)topelm=links;
    else{
      topelm=locate_prev(count-1);
      if(topelm)topelm=topelm->next;
    }
  } */

  lPair *last_link;
public:
  // constructors
  slList(){
    cmp_func=&default_cmp;
    last_link=NULL;
    //topelm=NULL;
  }

  slList(cmpf func){
    cmp_func=func;
    last_link=NULL;
    //topelm=NULL;
  }

  slList(const slList& other):lList<Item>(other){
    cmp_func=other.cmp_func;
    last_link=other.last_link;
    //locate_topelm();
  }

  slList& operator=(const slList& other){
    if(this!=&other){
      Remove_all();
      equate(other);
      cmp_func=other.cmp_func;
      last_link=other.last_link;
    }
    return *this;
  }

  virtual basicList<Item>* Clone() const{
    return new slList(*this);
  }

  // setting new sort function  (and sorting accordingly)
  virtual void SetFunc(cmpf func){
    cmp_func=func;
    SetList(this);
  }



  // sorted Append
  // inserts na items
  virtual long AppendP(const Item &a, long na=1, void **appelm=NULL);

  // gets the last insertion point
  virtual lPair *GetLastIns() const{
    return last_link;
  }

  // gets the element at last insertion position
  virtual Item *GetLastInsElm() const{
    lPair *pos=GetLastIns();
    if(pos)return &(pos->it);
    return NULL;
  }


  // sorted Iterate
  virtual void Iterate(Item (*func)(long,Item)){
    slList *copy=new slList(*this);
    Remove_all();
    long curp=copy->Rewind();
    while(!copy->Done()){
      Item a=func(copy->GetCurPos(),copy->GetCur());
      Append(a);
      copy->Next();
    }
    SetCurPos(curp);
    delete copy;
  }

  // attention: sorts after setting!!!
  virtual void SetCur(Item a){
    long curp=GetCurPos();
    Remove(curp);
    Append(a);
    SetCurPos(curp);
  };

  /*
  virtual Item &GetBottom() const {
    if(!links)fatal_error("Getting Bottom of void list!\n");
    else return links->it;
  }

  virtual Item &GetTop() const {
    if(!topelm)fatal_error("Getting Top of void list!\n");
    else return topelm->it;
  } */


};


template <class Item> Item *del_pointer(long i,Item *ptr){
  delete ptr;
  return NULL;
};


// ---------------- lList-------------------------------------


template <class Item> 
typename lList<Item>::lPair* lList<Item>::locate_prev(long index, lPair *first, long start_i) const{
  if(!first)first=links;
  lPair *elm, *prev=NULL;
  if(!first || index>=count)return NULL;
  else elm=first;
  long i=start_i;
  while(elm->next){
    if(i==index)return prev;
    prev=elm;
    elm=elm->next;
    i++;
  }
  return prev;
};


/*
template <class Item> lList<Item>::lList(const lList& other):basicList<Item>(other){
  equate(other);
}

template <class Item> lList& lList<Item>::operator=(const lList& other){
  Remove_all();
  equate(other);
} */

template <class Item> void lList<Item>::equate(const lList& other){
  count=other.count;
  if(!other.links){
    last=links=NULL;
    curelm=NULL;
  }
  else{
    // creating structure
    links = new lPair();
    if(!links)fatal_error("lList: MAE copying list!\n");
    //rewind_it();

    lPair *elmo=other.links, *elm=links;
    elm->it=elmo->it;
    while(elmo->next){
      elmo=elmo->next;
      elm->next=new lPair();
      if(!elm->next)fatal_error("lList: MAE copying list!\n");
      elm=elm->next;
      elm->it=elmo->it;
    }
    elm->next=NULL;
    last=elm;
    Rewind();
    SetCurPos(other.curpos);
  }
}


template <class Item> Item& lList<Item>::Get(long index) const {
  if(index<=0){
    if(!links)fatal_error("lList.Get: trying to Get from void list\n");
    else return links->it;
  }
  lPair *prev=locate_prev(index);
  if(!prev){
    fatal_error("lList.Get: index %ld out of range [0,%ld]\n",index,count);
  }
  return (prev->next)->it;
};


// returns index of the last appended element
template <class Item> long lList<Item>::AppendP(const Item& item, long na, void **appelm){
  if(na<1)fatal_error("slList: can't append %d elements!\n",na);
  lPair *elm;
  if(!links){
    links= new lPair;
    elm=links;
  }
  else{
    /*elm=links;
    while(elm->next)elm=elm->next;*/
    elm=last;
    elm->next= new lPair;
    elm= elm->next;
  }

  if(appelm)*appelm=(void *)elm;
  long i=na-1;
  do{
    if(!elm)fatal_error("lList: MAE appending element!\n");
    elm->it=item;
    count++;

    if(i==0){
      elm->next=NULL;
      last=elm;
      break;
    }
    elm->next= new lPair;
    elm=elm->next;
    i--;
  }while(1);

  return  count;
};


// returns new count
template <class Item> long lList<Item>::Remove(long index, long na){
  lPair *delptr, *elm, *prev=NULL;
  if(na<0)na=count;

  int rem_last=0;
  if(index+na>=count)rem_last=1; // last being removed

  if(!index){
   if(!links)fatal_error("lList: trying to Remove from void list!\n");
   elm=links;
   //links=elm->next;
  }
  else{
    prev=locate_prev(index);
    if(!prev){
      fatal_error("lList: Removing Item# %ld out of range [0,%ld]\n",index,count);
    }
    elm=prev->next;
  }


  long i=0;
  while(elm){
    if(i>=na)break;
    delptr=elm;
    elm=elm->next;
    delete delptr;
    count--;
    i++;
  }
  if(prev)prev->next=elm;
  else links=elm;

  if(rem_last){
    if(prev)last=prev;
    else last=links;
  }
  if(curpos>=index)rewind_it(); //if(curpos>=index)rewind_it();
  //if(curpos>count)curpos=count;
  return count;
};

// ---------------- lList-------------------------------------


// ---------------- slList-------------------------------------



template <class Item> long slList<Item>::AppendP(const Item &a, long na, void **appelm){
  if(na<1)fatal_error("slList: can't append %d elements!\n",na);
  long curp=this->Rewind();

  int i,res;
  lPair *prev=NULL, *telme,*telm0;
  // creating chunk
  for(i=0;i<na;i++){
    telme=  new lPair; // create_lPair();  //
    if(!telme)fatal_error("slList: MAE inserting element!\n");
    if(i==0){
      telm0=telme;
      if(appelm)*appelm=(void *)telm0;
    }
    telme->it=a;
    if(prev)prev->next=telme;
    prev=telme;
  }
  telme->next=NULL;

  long ind=-1;
  prev=NULL;
  while(!this->Done()){
    res=cmp_func(&(this->curelm->it),&a);
    if(res>0){
      telme->next=this->curelm;
      break;
    }
    prev=this->curelm;
    ind=this->Next();
  }
  if(res<=0){ // appended as last
    this->last=telme;
  }

  if(prev)prev->next=telm0;
  else this->links=telm0;
  this->count+=na;

  this->Rewind();
  SetCurPos(curp);
  last_link=prev;
  return ind+1;
}

// ---------------- slList-------------------------------------

# endif
