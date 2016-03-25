/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005-2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   This source code is Free Software and distributed under the terms of wxWidgets license (www.wxwidgets.org) 
 *
 *   $Revision: 1.8 $
 *   $Date: 2013/05/24 17:53:14 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/seqpack.h,v 1.8 2013/05/24 17:53:14 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/seqpack.h,v $
$Revision: 1.8 $
$Author: valuev $
$Date: 2013/05/24 17:53:14 $
*/
/*s****************************************************************************
 * $Log: seqpack.h,v $
 * Revision 1.8  2013/05/24 17:53:14  valuev
 * sync with kintech svn
 *
 * Revision 1.7  2012/09/07 14:26:06  valuev
 * removed using namespace from h files
 *
 * Revision 1.6  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.38  2010/10/07 10:30:22  lesha
 * restorer is commented
 *
 * Revision 1.37  2010/09/29 10:40:49  lesha
 * comments
 *
 * Revision 1.36  2010/04/17 19:32:39  valuev
 * fixed group_pack
 *
*******************************************************************************/


#ifndef SEQPACK_H
#define SEQPACK_H

/* \en @file seqpack.h
    @brief Classes for packing iterator sequences. Iterator increments are packed in vector<int>.
    
    Классы для упаковки данных. 
    Они отличаются друг от друга типом данных, для упаковки которых они предназначены, 
    способом их хранения, однако имеют общий интерфейс для работы.
    Начало записи осуществляется с помощью функции start_record.
    Далее последовательно вызывается функция record, аргументом которой является пакуемое значение. 
    По окончании записи вызывается end_record. Производить дозапись после вызова end_record нельзя.
    Однако можно начать запись сначала, предварительно вызвав функцию clear или start_record.

    К данным предоставляется последовательный доступ с помощью итераторов, 
    интерфейс работы с которыми похож на интерфейс работы с итераторами STL.
    Функции begin и end возвращают первый и следующий за последним итератор.
    У итераторов имеются операторы ++ (переход к следующему итератору), * (взятие значения),
    != (сравнение со следующим за последним итератором). Иллюстрацией служит следующий код, 
    в котором присходит перебор всех запакованных значений объекта packer класса packer_t:
    for(packer_t::iterator it=packer.begin(), e=packer.end(); it!=e; ++it)
      packer_t::value_t val=*it;

    В данный момент, похоже, у функций clear и start_record одинаковая реализация, и они дублируют друг друга...
**/

#include <vector>
#include <list>
#include <functional>
#include <map>
#include <algorithm>
#include <string>

#include "refobj.h"
#include "utiltl.h"



/// packer iterator group categories
/// no groups defined
template <class packer_it>
struct nogroup_it{
  typedef packer_it group_it;
  static group_it get_group_begin(const packer_it &it){
    return it;
  }
  static int get_group_count(const packer_it &it){
    return 1;
  }
};

template< class packer_it>
struct packit_traits{
  typedef typename packer_it::gcategory::group_it  group_it;
  typedef typename packer_it::gcategory gcategory;
};

/// default for all pointers
template< class T >
struct packit_traits< T* >{
  typedef T*         group_it;
  typedef nogroup_it<T *> gcategory;
};

/// extracts the next value at a given position 
/// pos must be less than or equal ind.size()-2
/// @returns the value or -1 if end of sequence or the postion is incorrect
/// modifies pos (the next reading position) and incrcount
int get_next_record(const vector<int> &ind, int &pos, int &incrcount, int single=0);

/// specifies the next value to be recorded
/// only nonnegative values are possible
/// @retuns    0 if the value was not packed (just stored)
///            1 packed with single increment
///            2 packed with other increment
///            3 packed as continued sequence
///           -1 tried to store after end-of sequence 
int put_next_record(vector<int> &ind, int cur, int single=0);

/// puts the end of record indicator 
int put_record_end(vector<int> &ind, int single=0);

/// class for packing integers
class int_pack{
  vector<int> ind;
  int single;
  size_t count;
public:
  typedef int value_t;
  typedef int record_t;

  int_pack(int sel_incr=0):single(sel_incr),count(0){
  }
  class iterator{
    friend class int_pack;
    int pos;
    int incrcount;
    int cur;
    const int_pack *parent;
    iterator(const int_pack *sparent, int end=0): pos(0),incrcount(0),parent(sparent) {
      if(!end){
        cur=get_next_record(parent->ind,pos,incrcount,parent->single);
        //cur=parent->ind[0];
      }
      else{
        pos=(int)parent->ind.size()-1;
        //while(pos>=0 && parent->ind[pos]<0)pos--;
        //pos++;
        cur=-1;
      }
    }
  public:
    typedef nogroup_it<iterator> gcategory;
    iterator(const iterator &other):pos(other.pos),incrcount(other.incrcount),cur(other.cur), parent(other.parent){
    }
    iterator() {} // for iterators array
    int operator*() const {
      return cur;
    }
    int operator++(){ // prefix
      return (cur=get_next_record(parent->ind,pos,incrcount,parent->single));
    }
    int operator++(int){ // postfix
      int tmp=cur;
      ++*this;
      return tmp;
    }
    /// checks for the end only
    bool operator!=(const iterator &other) const {
      if(pos!=other.pos)return true;
      else return false;
    }
  };
  int start_record(){
    ind.clear();
    count=0;
    return 1;
  }
  
 
  int next_record(int i){
    count++;
    return put_next_record(ind,i,single);
  }


   /// put a group of n elements starting with beg
  template <class int_it>
  int next_group(int_it beg, int n=1){
    int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    return put_record_end(ind,single);
  }

  iterator begin() const {
    if(ind.size()<3)return end();
    return iterator(this,0);
  }
  iterator end() const {
    return iterator(this,1);
  }

  /// returns size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof(int)*ind.size();
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(int)*size();
  }

  size_t size() const {
    return count;
  }

  void clear(){
    count=0;
    ind.clear();
  }
};


/// reads the input to int_pack from the string in the form of a list: 1-10,3,4,3-65
/// the integers must be positive inside a range [imin,imax]
/// the ranges num1-num2 are by default with increment 1 (-1)
/// delimiters (,-) are specified in the delims string
/// special value "all" means the whole range [imin,imax]
/// @return:  -1 wrong arguments,
///           -2 read (format) error
///           -3 outside the range
///            1 OK
int add_ranges(int_pack &packer, const string &spec, int imin=0, int imax=1000, bool finalize=true,const string &delims=",-");


/// translates a selected set of iterator values into an iterator
template <class iter_t>
class iter_pack{
  int_pack ipack;
  iter_t beg;
public:
  typedef iter_t value_t;
  typedef int record_t;
  typedef int_pack::iterator inc_iterator;
  
  iter_pack(iter_t begin, int sel_incr=0):beg(begin), ipack(sel_incr){}
  
  iter_pack(int sel_incr=0):ipack(sel_incr){}
  
  
  class iterator{
    friend class iter_pack;
  
    iter_t curit;
    int_pack::iterator incr;
    int_pack::iterator iend;

    iterator(const iter_pack *parent, int end=0):curit(parent->beg),incr(end ? parent->ipack.end() : parent->ipack.begin()), iend(parent->ipack.end()){
      if(!end && incr!=iend)curit+=*incr; // proceed to the first increment
    }
  public:
    typedef nogroup_it<iterator> gcategory;
    //typedef typename iter_t::value_t value_t;
    iterator(const iterator &other):curit(other.curit),incr(other.incr),iend(other.iend){}
    iter_t operator*() const {
      return curit;
    }
    iterator &operator++(){ // prefix
      ++incr;
      if(incr!=iend)curit+=*incr;
      return *this;
    }
    iterator operator++(int){ // postfix
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    /// checks for the end only
    bool operator!=(const iterator &other) const {
      if(incr!=other.incr)return true;
      else return false;
    }

  };

  int start_record(){
    return ipack.start_record();
  }
  int next_record(int i){
    return ipack.next_record(i);
  }

  /// put a group of n elements starting with beg
  template <class int_it>
  int next_group(int_it beg, int n=1){
    return ipack.next_group(beg,n);
  }

  int end_record(){
    return ipack.end_record();
  }

  void set_start_it(iter_t sbeg){
    beg=sbeg;
  }

  iterator begin() const {
    return iterator(this);
  }
  iterator end() const {
    return iterator(this,1);
  }

  inc_iterator inc_begin() const {
    return ipack.begin();
  }

  inc_iterator inc_end() const {
    return ipack.end();
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return ipack.packed_size();
  }

  size_t size() const {
    return ipack.size();
  }
  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }

  void clear(){
    ipack.clear();
  }

};


template<class pair_t>
struct pair_sort_t {
  bool operator()(const pair_t& a, const pair_t& b) const {
    return a.first<b.first;
  }
};

/// reorders two packs of the same size in memory with < operator, the first is pack the leading, 
/// the second follows the order of the first
template<class pack_t1, class pack_t2>
int mem_reorder(pack_t1 &p1, pack_t2 &p2){
  size_t n=p1.size(), i;
  if(n!=p2.size())return -1;
  typename pack_t1::iterator it1=p1.begin(), e1=p1.end();
  typename pack_t2::iterator it2=p2.begin(), e2=p2.end();
  typedef typename pack_t1::value_t value_t1;
  typedef typename pack_t2::value_t value_t2;
  typedef pair<value_t1, value_t2> sort_t;

  pair_sort_t<sort_t> pred;
  
  // copying
  sort_t *v= new sort_t[n];
  for(i=0;i<n;i++){
    v[i]=sort_t(*it1,*it2);
    ++it1;
    ++it2;
  }
  // sorting
  sort(v,v+n,pred);
  // copying the sorted records
  p1.clear();
  p1.start_record();
  p2.clear();
  p2.start_record();
  for(i=0;i<n;i++){
    p1.next_record(v[i].first);
    p2.next_record(v[i].second);
  }
  p1.end_record();
  p2.end_record();
  delete [] v;
  return 1;
}


/// function to test integer sequence packing
int test_seqpack(int len,  double seqprob, int seqlenav, int spec_val, double specprob);

template <class comp_pr>
struct subgroup_test{
  /// tests whether the begining of sequence pointed by beg1 coincides with the sequence 
  /// of n elements pointed by beg2
  /// @return the number of non-coinciding elements (returns 0 if beg2 is a subsequence) 
  template <class g1_it, class g2_it >
  int operator()(g1_it beg1, g1_it end1, g2_it beg2, int n, const comp_pr &pr) const {
    while(n>0 && beg1!=end1){
      if(!pr(*beg1,*beg2))break;
      ++beg1;
      ++beg2;
      n--;
    }
    return n;
  }
};


/// container with duplicates control using comparison prediacate of type comp_pr
/// which must have bool opertor()(const T&,const T&) returning true for duplicated items
/// supports vectors or lists as container types
template<class T, class cont_t=vector<T>, class comp_pr=equal_to<T> >
class container_set: public cont_t{
protected:
  comp_pr pr;
public:
  typedef typename cont_t::iterator iterator;
  typedef comp_pr key_compare;
  /// returns the iterator of existing element if it is present (checks by operator==) and false
  //  otherwise 
  /// appends element to the end of the vector and returns last element iterator and  true  
  pair<iterator, bool> insert(const T& elm){
    iterator it=cont_t::begin(), e=cont_t::end(); 
    for(;it!=e;++it){
      if(pr(*it,elm))return make_pair(it,false);
    }
    push_back(elm);
    return make_pair(e,true); // is that OK?
  }

  /// returns the number of elements between begin and the inserted element
  size_t insert_ind(iterator it,const T& elm){
    iterator e=cont_t::end(); 
    size_t i=0;
    for(;it!=e;++it){
      if(pr(*it,elm))return i;
      i++;
    }
    push_back(elm);
    return i; 
  }


  /// the same with default starting
  size_t insert_ind(const T& elm){
    return insert_ind(cont_t::begin(),elm);
  }

  template<class inp_it, class group_pr>
  size_t insert_group(iterator it, inp_it beg, int n, const group_pr & gpr){
    iterator e=cont_t::end(); 
    size_t i=0;
    for(;it!=e;++it){
      if(gpr(it,e,beg,n,pr)==0){ // found subsequence which is exactly n elements long
        return i;
      }
      i++;
    }
    // not found: inserting
    while(n>0){
      push_back(*beg);
      ++beg;
      n--;
    }
    return i;
  }

  template<class inp_it, class group_pr>
  size_t insert_group(inp_it beg, int n, const group_pr & gpr){
    return insert_group(cont_t::begin(), beg, n, gpr);
  }

  void clear_index(){}
};


template<class T, class comp_pr=equal_to<T> > 
class vector_set: public container_set<T,vector<T>, comp_pr >{
public:
  typedef typename container_set<T,vector<T>, comp_pr >::iterator iterator;
  typedef typename container_set<T,vector<T>, comp_pr >::key_compare key_compare;
};



template<class T, class comp_pr=equal_to<T> > 
class list_set: public container_set<T,list<T>, comp_pr >{
public:
  typedef typename container_set<T,list<T>, comp_pr >::iterator iterator;
  typedef typename container_set<T,list<T>, comp_pr >::key_compare key_compare;
};

template<class T, class less_pr>
struct equal_from_less_pr: public binary_function<T,T,bool>{
  less_pr pr;
  bool operator()(const T& r1, const T& r2) const {
    return (!pr(r1,r2))&&(!pr(r2,r1));
  }
};

template<class T, class less_pr=less<T> > 
class vector_indexed: public container_set<T,vector<T>, equal_from_less_pr<T,less_pr> >{
public:
  typedef container_set<T,vector<T>, equal_from_less_pr<T,less_pr> > base_t;
  //typedef container_set<T,vector<T>, equal_to<T> /*equal_from_less_pr<T,less_pr>*/ > base_t;
  typedef typename base_t::iterator iterator;
  typedef typename base_t::key_compare key_compare;
protected:
  typedef multimap<T, size_t, less_pr> map_t;
  map_t imap;
  typedef typename map_t::iterator map_it;
  typedef typename map_t::value_type map_vt;

public:

  template<class inp_it, class group_pr>
  size_t insert_group(typename base_t::iterator it, inp_it beg, int n, const group_pr & gpr){
    size_t i=0;
    map_it mit=imap.find(*beg);
    // check the group in any case
    iterator e=this->end(); 
    while(it!=e){
      if(mit!=imap.end() && mit->first==*beg){
        it=base_t::begin()+mit->second;
        i=mit->second;
        ++mit;
      }
      if(gpr(it,e,beg,n,this->pr)==0){ // found subsequence which is exactly n elements long
        return i;
      }
      i++;
      ++it;
    }
    int i0=(int)i;
    // not found: inserting
    while(n>0){
      push_back(*beg);
      imap.insert(map_vt(*beg,i0));
      ++beg;
      n--;
      i0++;
    }
    return i;
  }

  template<class inp_it, class group_pr>
  size_t insert_group(inp_it beg, int n, const group_pr & gpr){
    return insert_group(this->begin(), beg, n, gpr);
  }


  /// insert with it starting
  size_t insert_ind(typename base_t::iterator it, const T& elm){
    map_it mit=imap.find(elm);
    if(mit!=imap.end())return mit->second;
    // inserting
    imap.insert(map_vt(elm,base_t::size()));
    push_back(elm);
    return base_t::size()-1;
  }
  /// the same with default starting
  size_t insert_ind(const T& elm){
    return insert_ind(base_t::begin(),elm);
  }

  void clear_index(){
    imap.clear();
  }

  void clear(){
    clear_index();
    base_t::clear();
  }

};



/// class for unpacked data stored in non-packing container of cont_t 
/// that emulates data_pack
template<class T, class cont_t=vector<T> >
class data_unpacked{
public:
  typedef typename cont_t::iterator iterator;
  typedef typename cont_t::iterator vset_it;
  typedef T value_t;
  typedef T record_t;

protected: 
  cont_t vset;
public:
  
  data_unpacked(){}

  int start_record(){
    vset.clear();
    return 1;
  }
  int next_record(const T& value){
    vset.push_back(value);
    return 1;
  }

   /// put a group of n elements starting with beg
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }


  int end_record(){
    return 1;
  }
  
  iterator begin(){
    return vset.begin();
  }
  
  iterator end(){
    return vset.end();
  }

  size_t size() const {
    return vset.size();
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return data_size();
  }
  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }

  void clear(){
    vset.clear();
  }
};
#if 0
template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr> >
class data_vset{
public:
  typedef typename set_t::iterator vset_it;
  typedef T value_t;
  typedef T record_t;

  //friend class data_vset<T, set_t>::iterator;
protected: // ! if protected thereis error in group_pack 
//public:
  set_t vset;
public:
  vector<int> ipack;
  class iterator{
  protected:
    friend class data_vset;
    vector<int>::iterator it;
    data_vset *parent;
    iterator(vector<int>::iterator sit, data_vset *sparent): it(sit), parent(sparent){}
  public:  
    typedef nogroup_it gcategory;
    /// copy constructor
    iterator(const iterator &other):it(other.it), parent(other.parent){}
    iterator() {} // for iterators array
    
    /// WARNING: recording by the given reference will  effectively
    ///          change ALL data entries that are equal to this one 
    T& operator*() const {
      return parent->vset[*it];
    }
    // prefix increment
    iterator& operator++(){
      ++it;
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++it;
      return tmp;
    } 

    bool operator!=(const iterator &other){
      return (it!=other.it);
    } 

  };
  data_vset(): vset(), ipack(0) {}

  int start_record(){
    vset.clear();
    ipack.clear();
    return 1;
  }
  int next_record(const T& value){
    size_t res=vset.insert_ind(value); //(size_t)(vset.insert(value).first-vset.begin());
    ipack.push_back((int)res);
    return 1;
  }

  /// put a group of n elements starting with beg
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    vset.clear_index();
    return 1;
  }
  
  iterator begin() {
    return iterator(ipack.begin(),this);
  }
  
  iterator end()  {
    return iterator(ipack.end(),this);
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof(value_t)*vset.size();
  }

  size_t size() const {
    return ipack.size();
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }
  void clear(){
    vset.clear();
    ipack.clear();
  }
};
#endif
/// Class for packing data using packing container cont_t.
/// cont_t can be container_set<T, container_t>, or just set<T,sorter_t>: NOW WORKS WITH vector_set or list_set only.
/// The comparison predicate is of type comp_pr
/// which must have bool opertor()(const T&,const T&) returning true for duplicated items
template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr> >
class data_pack{
public:
  typedef typename set_t::iterator vset_it;
  typedef T value_t;
  typedef T record_t;

  //friend class data_pack<T, set_t>::iterator;
protected: // ! if protected thereis error in group_pack 
//public:
  set_t vset;
public:
  int_pack ipack;
  class iterator{
  protected:
    friend class data_pack;
    int_pack::iterator it;
    data_pack *parent;
    iterator(int_pack::iterator sit, data_pack *sparent): it(sit), parent(sparent){}
  public:  
    typedef nogroup_it<iterator> gcategory;
    /// copy constructor
    iterator(const iterator &other):it(other.it), parent(other.parent){}
    iterator() {} // for iterators array
    
    /// WARNING: recording by the given reference will  effectively
    ///          change ALL data entries that are equal to this one 
    T& operator*() const {
      return parent->vset[*it];
    }
    // prefix increment
    iterator& operator++(){
      ++it;
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++it;
      return tmp;
    } 

    bool operator!=(const iterator &other){
      return (it!=other.it);
    } 

  };
  data_pack(): vset(), ipack(0) {}

  int start_record(){
    vset.clear();
    return ipack.start_record();
  }
  int next_record(const T& value){
    size_t res=vset.insert_ind(value); //(size_t)(vset.insert(value).first-vset.begin());
    return ipack.next_record((int)res);
  }

  /// put a group of n elements starting with beg
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }

  int end_record(){
    vset.clear_index();
    return ipack.end_record();
  }
  
  iterator begin() {
    return iterator(ipack.begin(),this);
  }
  
  iterator end()  {
    return iterator(ipack.end(),this);
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof(value_t)*vset.size()+ipack.packed_size();
  }

  size_t size() const {
    return ipack.size();
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size();
  }
  void clear(){
    vset.clear();
    ipack.clear();
  }
};

/// Class for packing pairs of data of two different types.
/// The types are packed by different packers but recovered simultaneously by an iterator
/// pair_tt is the returned value type, it must have at least one of the following constructors
///          pair_tt(value_t1 &a,value_t2 &b) 
///          pair_tt(const value_t1 &a,const value_t2 &b)
///          pair_tt(value_t1 a,value_t2 b) 
/// which may construct a type passing arguments either by value or by reference. 
/// WARNING: recording by the returned references is possible for some packers, but may lead to changes 
///          in other packed elements,
///          the policy of recording is determined by packers
/// By default the returned value type of pair_pack is pair<value_t1,value_t2> returned by value (no recording possible)
/// see refpair_pack()
template<class T1, class T2, class packer_t1=data_pack<T1>, class packer_t2=data_pack<T2>, 
         class pair_tt=pair<typename packer_t1::value_t,typename packer_t2::value_t > >
class pair_pack{
public: 
  packer_t1 *pack1;
  packer_t2 *pack2;
 

  typedef typename packer_t1::record_t record_t1;
  typedef typename packer_t2::record_t record_t2;
  typedef pair_tt value_t;
  typedef pair<record_t1,record_t2> record_t;
  
  //friend class pair_pack<T1,T2,packer_t1,packer_t2>::iterator;
  
  /// constructor with possible packer construction specification
  pair_pack(packer_t1 *sp1=new packer_t1(), packer_t2 *sp2 = new packer_t2()):pack1(sp1),pack2(sp2){}

  class iterator{
    friend class pair_pack;
    typename packer_t1::iterator it1;
    typename packer_t2::iterator it2;
    iterator(typename packer_t1::iterator sit1, typename packer_t2::iterator sit2): it1(sit1),it2(sit2){}
  public:  
    /// the class gcategory must have a function
    /// static int get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;

    static int get_group_count(const iterator &it){ 
      int n1=get_group_count(it.it1);
      int n2=get_group_count(it.it2);
      return (n1>n2 ? n1 : n2);
    }

    /// copy constructor
    iterator(const iterator &other):it1(other.it1),it2(other.it2){}
    iterator() {} // for iterators array
    
    /// WARNING: recording by the given references will  effectively
    ///          change ALL data entries that are equal to the referenced ones 
    value_t operator*() const {
      return value_t(*it1,*it2);
    }
    /// prefix increment
    iterator& operator++(){
      ++it1;
      ++it2;
      return *this;
    } 
    /// postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 
    bool operator!=(const iterator &other){
      return (it1!=other.it1);
    }
  };
  
  int start_record(){
    int res=pack1->start_record();
    if(res<0)return res;
    res=pack2->start_record();
    return res;
  }
  
  int next_record(const record_t1 &v1, const record_t2 &v2){
    int res=pack1->next_record(v1);
    if(res<0)return res;
    res=pack2->next_record(v2);
    return res;
  } 
  
  int next_record(const record_t &p){
    return next_record(p.first,p.second);
  }

  /// put a group of n elements starting with beg
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    int res=pack1->next_group(first_it<rec_it>(beg), n);
    if(res<0)return res;
    res=pack2->next_group(second_it<rec_it>(beg), n);
    return res;
    /*int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;*/
  }
  
  int end_record(){
    int res=pack1->end_record();
    if(res<0)return res;
    res=pack2->end_record();
    return res;
  }

  iterator begin(){
    return iterator(pack1->begin(),pack2->begin());
  }
  
  iterator end(){
    return iterator(pack1->end(), pack2->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return pack1->packed_size()+pack2->packed_size();
  }

  size_t size() const {
    return pack1->size(); // must be the same for both
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return pack1->data_size()+pack2->data_size();
  }

  void clear(){
    pack1->clear();
    pack2->clear();
  }

  ~pair_pack(){
    delete pack1;
    delete pack2;
  }
};

/// both types stored by reference
template <class T1, class T2>
struct refpair{
  T1& first;
  T2& second;
  typedef T1 first_type;
  typedef T2 second_type;
  refpair(T1 &a, T2 &b):first(a),second(b){}

  operator pair<T1,T2>(){
    return make_pair(first,second);
  }

  T1 &r1() const { return first; }
  T2 &r2() const { return second; }
  refpair(const refpair<T1,T2> &other): first(other.first), second(other.second) {}
};

/// the first type is stored by value, the second by reference
template <class T1, class T2>
struct vrefpair{
  T1  first;
  T2& second;
  typedef T1 first_type;
  typedef T2 second_type;
  vrefpair(const T1 &a, T2 &b):first(a),second(b){}
  T2 &r2() const { return second; }
  vrefpair(const vrefpair<T1,T2> &other): first(other.first), second(other.second) {}
};

/// the first type is stored by reference, the second by value
template <class T1, class T2>
struct refvpair{
  T1&  first;
  T2 second;
  typedef T1 first_type;
  typedef T2 second_type;
  refvpair(T1 &a, const T2 &b):first(a),second(b){}
  T1 &r1() const { return first; }
  refvpair(const refvpair<T1,T2> &other): first(other.first), second(other.second) {}
};


template<class T1, class T2, class packer_t1=data_pack<T1>, class packer_t2=data_pack<T2> > 
class refpair_pack: public pair_pack<T1,T2,packer_t1,packer_t2, refpair<typename packer_t1::value_t,typename packer_t2::value_t> >{
public:
  typedef pair_pack<T1,T2,packer_t1,packer_t2, refpair<typename packer_t1::value_t,typename packer_t2::value_t> > base_t;
  typedef typename base_t::record_t1 record_t1;
  typedef typename base_t::record_t2 record_t2;
  typedef typename base_t::value_t value_t;
  typedef typename base_t::record_t record_t;
  typedef typename base_t::iterator iterator;

  /// constructor with possible packer construction specification
  refpair_pack(packer_t1 *sp1=new packer_t1(), packer_t2 *sp2 = new packer_t2()):base_t(sp1,sp2){}

};






/// Class for packing pairs (iterator,data) of two different types.
/// The iterator is icremented (not packed) and the data is packed by packer_t.
template<class iter_t, class T, class packer_t=data_pack<T> >
class seq_data_pack{
protected: 
  iter_t *beg;
  packer_t *pack;
public:
  typedef typename packer_t::record_t record_t;
  typedef pair<typename iter_t::value_t,typename packer_t::value_t> value_t;
  
  
  /// constructor with possible packer construction specification
  seq_data_pack(iter_t *sbeg=NULL, packer_t *sp=new packer_t()):pack(sp),beg(NULL){
    if(sbeg)set_start_it(*sbeg);
  }

  class iterator{
    friend class seq_data_pack;
    iter_t it1;
    typename packer_t::iterator it2;
    iterator(iter_t sit1, typename packer_t::iterator sit2): it1(sit1),it2(sit2){}
  public:  
    /// the class gcategory must have a function
    /// static int get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;

    static int get_group_count(const iterator &it){ 
      int n2=get_group_count(it.it2);
      return n2;
    }

    /// copy constructor
    iterator(const iterator &other):it1(other.it1),it2(other.it2){}
    
    /// WARNING: recording by the given references will  effectively
    ///          change ALL data entries that are equal to the referenced ones 
    value_t operator*() const {
      return value_t(*it1,*it2);
    }
    /// prefix increment
    iterator& operator++(){
      ++it1;
      ++it2;
      return *this;
    } 
    /// postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 
    bool operator!=(const iterator &other){
      return (it2!=other.it2);
    }
  };

  int set_start_it(iter_t sbeg){
    if(beg)delete beg;
    beg = new iter_t(sbeg);
    return 1;
  }
  
  int start_record(){
    return pack->start_record();
  }
  
  int next_record(const typename packer_t::record_t &v){
    return pack->next_record(v);
  } 

  /// put a group of n elements starting with beg
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    int i;
    for(i=0;i<n;++i,++beg){
      if(next_record(*beg)<0)break;
    }
    return i;
  }
   
  int end_record(){
    return pack->end_record();
  }

  iterator begin() const {
    return iterator(*beg,pack->begin());
  }

  iterator end() const {
    return iterator(*beg,pack->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return pack->packed_size()+sizeof(beg);
  }

  size_t size() const {
    return pack->size(); // must be the same for both
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(T)*pack->size()+sizeof(beg);
  }

  void clear(){
    pack->clear();
  }

  ~seq_data_pack(){
    delete pack;
    if(beg)delete beg;
  }
};


/// Class for packing iterator together with some packed data.
/// The iterator and data are packed by different packers but recovered simultaneously by an iterator.
/// By default the return type is pair<it,value>, recording by reference is not possible
template<class packed_it, class T,  class packer_t=data_pack<T>, class pair_tt= pair<packed_it, typename packer_t::value_t> >
class iter_data_pack : public pair_pack<packed_it, T, iter_pack<packed_it>,packer_t, pair_tt >{
public:
  typedef pair_pack<packed_it, T, iter_pack<packed_it>,packer_t, pair_tt  > base_t;
  typedef typename base_t::record_t record_t;
  typedef typename base_t::value_t value_t;
  
  
  typedef typename base_t::iterator iterator;
  
  /// constructor with possible packer construction specification
  iter_data_pack(packed_it beg, packer_t *sp=new packer_t()): base_t(new iter_pack<packed_it>(beg),sp){}
};


/// returns the number elements left in the sequence grouped with a given iterator
/// the element pointed by given iterator is also counted, so the default is 1
template <class packer_it>
inline int get_group_count(const packer_it &it){
  typedef typename packit_traits<packer_it>::gcategory cat_t;
  return cat_t::get_group_count(it);
}

/// returns the begining of the group associated with current iteratortemplate <class packer_it>
template <class packer_it>
inline typename packit_traits<packer_it>::group_it get_group_begin(const packer_it &it){
  typedef typename packit_traits<packer_it>::gcategory cat_t;
  return cat_t::get_group_begin(it);
}


/// Class for packing groups of data using packing container set_t.
/// cont_t can be container_set<T, container_t>, or just set<T,sorter_t>: NOW WORKS WITH vector_set or list_set only.
/// The comparison predicate is of type comp_pr
/// which must have bool opertor()(const T&,const T&) returning true for duplicated items
template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr> , class group_pr=subgroup_test<typename set_t::key_compare> >
class group_pack: public data_pack<T,comp_pr,set_t>{
public:
  typedef data_pack<T,comp_pr,set_t> base_t;
  typedef typename base_t::vset_it vset_it;
  typedef T value_t;
  typedef T record_t;
  /// STL-like name
  typedef group_pr key_compare;

  //friend class data_pack<T, set_t>::iterator;
protected: 
  using  base_t::vset;
  using  base_t::ipack;

  mngptr<int_pack> npack;
  group_pr pr;
  int sz;
public:

  class iterator: public base_t::iterator {
    typedef typename base_t::iterator base_it;
    friend class group_pack;
    int_pack::iterator itn;
    int n;
//    using typename base_it::it;
    using base_it::it;
    iterator(int_pack::iterator siti, group_pack *sparent, int_pack::iterator sitn): base_it(siti,sparent),itn(sitn),n(0){
    }
  public:  
    /// the class gcategory must have a function
    /// static int get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;
    typedef vset_it group_it;

    static int get_group_count(const iterator &it){ 
      return *(it.itn)-it.n;
    }

    static group_it get_group_begin(const iterator &it){ 
/* old version that worked OK: 6 hours of MPI debugging to catch this bug 
      if (it.it!=it.parent->ipack.end()) {
        int sh=*(it.it);
        group_pack *ig=(group_pack *)it.parent;
        return (ig)->vset.begin()+sh+it.n;
      }
      return ((group_pack *)it.parent)->vset.end(); */
       
      int sh=*(it.it);
      if(sh>=0){ // may be negative if gropus with 0 elements exist at the end of the list and *it is called!
        group_pack *ig=(group_pack *)it.parent;
        return (ig)->vset.begin()+sh+it.n;
      }
      return ((group_pack *)it.parent)->vset.end();
    }

    /// copy constructor
    iterator(const iterator &other):base_it(other),itn(other.itn),n(other.n){}
    iterator() {} // for iterators array

    /// WARNING: recording by the given reference will  effectively
    ///          change ALL data entries that are equal to this one 
    T& operator*() const {
      return ((group_pack *)this->parent)->vset[*it+n];
    }
    // prefix increment
    iterator& operator++(){
      n++;
      if(n>=*itn){ // safe also after end of sequence, returns -1
        if(*itn!=0)
          ++it;
        ++itn;
        n=0;
      }
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    } 

    iterator& plus(){
      if(*itn)
        ++it;
      ++itn;
      n=0;
      return *this;
    }

    bool operator!=(const iterator &other) const {
//      return ((it!=other.it)&&(itn!=other.itn)&&(n!=other.n));
      return (it!=other.it);
    } 

  };
  
  group_pack():sz(0),npack(new int_pack(1),1){}

  int start_record(){
    sz=0;
    if(npack.managed()){ // if this is not dependent service packer
      npack->start_record();
    }
    return base_t::start_record();
  }

  int next_record(const T& value){
    return next_group(&value,1);
  }

  /// this allows to share auxiliary data with other packer
  /// THE PACKERS MUST BE SYNCRONIZED EXTERNALLY: next_group must be called in sync with parent packer
  /// USE WITH CARE
  template <class packer_t>
  int aux_depends_on(const packer_t &other){
    npack.reset(other.get_grouper(),0);
    return 1;
  }

  int_pack *get_grouper() const {
    return npack.ptr();
  }

  /// put a group of n elements starting from beg
  /// packing the group with 0 elements will lead to the following behaviour when unpacked:
  /// get_group_count(it) will return 0, *it will return the next recorded element (if it exists)
  /// then after ++it, *it will return the same element, but get_group_count(it) will return the next group count
  /// THE ELEMENTS POINTED BY BEG WILL NOT BE RECORDED WHEN n=0
  /// DO NOT USE *it if get_group_count(it) returns 0!
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    if(npack.managed()){ // if this is not dependent service packer
      if(npack->next_record(n)<0)return -1;
    }
    if(n==0)
      return 0;
    size_t res=vset.insert_group(beg,n,pr);
    ipack.next_record((int)res);
    sz+=n;
    return n;
  }

  int end_record(){
    if(npack.managed()){ // if this is not dependent service packer
      npack->end_record();
    }
    return base_t::end_record();
  }
  
  iterator begin() {
    return iterator(ipack.begin(),this, npack->begin());
  }
  
  iterator end()  {
    return iterator(ipack.end(),this, npack->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return sizeof base_t::packed_size()+(npack.managed() ? npack->packed_size() :0);
  }

  size_t size() const {
    return sz;
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return sizeof(value_t)*size()+(npack.managed() ? npack->data_size() :0);
//    return (npack.managed() ? npack->data_size() :0);
  }

  void clear(){
    npack->clear();
    sz=0;
    return base_t::clear();
  }

};


template<class T>
class group_unpacked: public data_unpacked<T>{
public:
  typedef data_unpacked<T> base_t;
  typedef typename base_t::vset_it vset_it;
  typedef T value_t;
  typedef T record_t;

protected: 
  using  base_t::vset;

  mngptr<data_unpacked<int> > npack;
  int sz;
public:

  class iterator: public base_t::iterator {
    typedef typename base_t::iterator base_it;
    friend class group_unpacked;
    data_unpacked<int>::iterator itn;
    int n;
    iterator(base_it siti, data_unpacked<int>::iterator sitn): base_it(siti),itn(sitn),n(0){}
  public:  
    /// the class gcategory must have a function
    /// static int get_group_count(const iterator &) which returns group count
    typedef iterator gcategory;
    typedef base_it group_it;
    static int get_group_count(const iterator &it){ 
      return *(it.itn)-it.n;
    }

    static group_it get_group_begin(const iterator &it){ 
      return base_it(it);
    }

    iterator() {} // for iterators array

    iterator& operator++(){
      n++;
      if(n>=*itn){ // safe also after end of sequence, returns -1
        if(*itn!=0)
          base_it::operator++();
        ++itn;
        n=0;
      }
      else
        base_it::operator++();
      return *this;
    } 
    // postfix increment
    iterator operator++(int){
      iterator tmp=*this;
      ++*this;
      return tmp;
    }
    iterator& plus(){
      base_it::operator+=(*itn-n);
      ++itn;
      return *this;
    }

    bool operator!=(const iterator &other) const {
      return itn!=other.itn;
    } 

  };
  
  group_unpacked():sz(0),npack(new data_unpacked<int>(),1){}

  int start_record(){
    sz=0;
    if(npack.managed()){ // if this is not dependent service packer
      npack->start_record();
    }
    return base_t::start_record();
  }

  int next_record(const T& value){
    return next_group(&value,1);
  }

  /// this allows to share auxiliary data with other packer
  /// THE PACKERS MUST BE SYNCRONIZED EXTERNALLY: next_group must be called in sync with parent packer
  /// USE WITH CARE
  template <class packer_t>
  int aux_depends_on(const packer_t &other){
    npack.reset(other.get_grouper(),0);
    return 1;
  }

  data_unpacked<int> *get_grouper() const {
    return npack.ptr();
  }

  /// put a group of n elements starting from beg
  /// packing the group with 0 elements will lead to the following behaviour when unpacked:
  /// get_group_count(it) will return 0, *it will return the next recorded element (if it exists)
  /// then after ++it, *it will return the same element, but get_group_count(it) will return the next group count
  /// THE ELEMENTS POINTED BY BEG WILL NOT BE RECORDED WHEN n=0
  /// DO NOT USE *it if get_group_count(it) returns 0!
  template <class rec_it>
  int next_group(rec_it beg, int n=1){
    if(npack.managed()){ // if this is not dependent service packer
      if(npack->next_record(n)<0)return -1;
    }
    if(n==0)return 0;
    base_t::next_group(beg,n);
    sz+=n;
    return n;
  }

  int end_record(){
    if(npack.managed()){ // if this is not dependent service packer
      npack->end_record();
    }
    return base_t::end_record();
  }
  
  iterator begin() {
    return iterator(base_t::begin(), npack->begin());
  }
  
  iterator end()  {
    return iterator(base_t::end(), npack->end());
  }

  /// returns the size in bytes ignoring vector overheads
  size_t packed_size() const {
    return data_size();
  }

  size_t size() const {
    return sz;
  }

  /// returns full unpacked sequence size in bytes 
  size_t data_size() const {
    return base_t::data_size()+(npack.managed() ? npack->data_size():0);
//    return (npack.managed() ? npack->data_size():0);
  }

  void clear(){
    npack->clear();
    sz=0;
    return base_t::clear();
  }
};

#endif
