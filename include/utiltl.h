# ifndef _UTILTL
# define _UTILTL

/*e @file utiltl.h
    @brief Useful templates (iterators, modifiers, etc)
  
**/

#include <iterator>
#include <complex>
#include <string.h>
#include "vector_3.h"
#include "logexc.h" 

template<class iter>
struct value_type_cast{
  typedef typename iter::value_type value_type;
};



template <class inp_it>
size_t string_scan_t(const char *str,char *format,char *delim, inp_it beg, inp_it end){
  char *buf=new char[strlen(str)+1];
  strcpy(buf,str);
  typename iterator_traits<inp_it>::value_type val;
  char *c1=strtok(buf,delim);
  size_t c=0;
  while(c1 && beg!=end){
    int res=sscanf(c1,format,&val);
    if(res>0){
      c++;
      if(beg!=end)
        *beg++=val;
      else
        break;
    }
    c1=strtok(NULL,delim);
  }
  delete [] buf;
  return c;
}


struct void_check_t{
  bool operator()(int i) const { (void)i; return true; }
};

/// aggregate iterator.
template<class it1, class it2=it1, class value_t= typename value_type_cast<it1>::value_type>
class aggregate_it{
  std::pair<it1,it1> p1;
  it2 i2;
public:
  aggregate_it(it1 beg1, it1 end1, it2 beg2):p1(beg1,end1),i2(beg2){}

  aggregate_it(const aggregate_it &other):p1(other.p1),i2(other.i2){}

  typedef value_t value_type;

  const value_type operator*() const {
    if(p1.first!=p1.second)return *p1.first;
    else return *i2;
  }
  
  aggregate_it& operator++(){ // prefix
    if(p1.first!=p1.second)p1.first++;
    else i2++;
    return *this;
  }

  aggregate_it operator++(int){ // postfix
    aggregate_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const aggregate_it &other) const{
    return (p1.first!=other.p1.first) || (i2!=other.i2);
  }
};

/// iterator removing some subsequence (beg2,end2) from the basic sequence (beg1,end1) if comparison predicate returns false
template<class iter1, class iter2=iter1, class pred_t=not_equal_to<iter1> >
class remove_if_it{
  iter1 beg1, end1;
  iter2 beg2, end2;
  pred_t Pr;
  
  // протаскивает beg1 сквозь последовательность (beg2,end2), если beg1 в нее попадает
  void check_beg1(){
    iter2 it2=beg2;
    for(;beg1!=end1 && it2!=end2;++it2){
      if(!Pr(beg1,it2))++beg1;
    }
  }

public:
  remove_if_it(iter1 sbeg1, iter1 send1, iter2 sbeg2, iter2 send2):beg1(sbeg1),beg2(sbeg2),end1(send1),end2(send2){
    check_beg1();
  }
  remove_if_it(const remove_if_it &other):beg1(other.beg1),beg2(other.beg2),end1(other.end1),end2(other.end2){
    check_beg1();
  }

  typedef typename value_type_cast<iter1>::value_type value_type;

  const value_type operator*() const {
    return *beg1;
  }
  
  remove_if_it& operator++(){ // prefix
    ++beg1;
    check_beg1();    
    return *this;
  }

  remove_if_it operator++(int){ // postfix
    remove_if_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator!=(const remove_if_it &other) const{
    return (beg1!=other.beg1);
  }
};

///\en  Skips some of the numbers, returned as false by checker
///     Checker MUST return true for out_of range indicies.
template<class iterator, class checker_t=void_check_t>
struct skiper_it{
  checker_t checker;
  iterator it;
  int i;
  
  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
	typedef typename iterator_traits<iterator>::value_type value_type;
	typedef typename iterator_traits<iterator>::difference_type difference_type;
	typedef difference_type distance_type;	// retained
	typedef typename iterator_traits<iterator>::pointer pointer;
	typedef typename iterator_traits<iterator>::reference reference;

  skiper_it(){}

  ///\en Use \a valid=flase to indicate invalid iterators, 
  ///    where no checks are performed (end of sequence, etc.)
  skiper_it(iterator it_, int i_=0, bool valid=true):it(it_),i(i_){
    if(valid){
      while(!checker(i)){
        ++i;
        ++it;
      }
    }
  }

  ///\en Use \a valid=flase to indicate invalid iterators, 
  ///    where no checks are performed (end of sequence, etc.)
  skiper_it(iterator it_, const checker_t &ch, int i_=0, bool valid=true):it(it_),i(i_),checker(ch){
    if(valid){
      while(!checker(i)){
        ++i;
        ++it;
      }
    }
  }
  
  
  value_type &operator*(){
    return *it;
  }
  skiper_it &operator++(){
    do{
      ++i;
      ++it;
    }while(!checker(i));
    return *this;
  }
  
  skiper_it operator++(int){ // postfix
    skiper_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const skiper_it &other) const {
    return (it==other.it);
  }
  bool operator!=(const skiper_it &other) const {
    return !(*this==other);
  }

};



/// looper: loops the sequence after nd iterations, counts the loop number (iter)
template<class iterator, class checker_t=void_check_t>
struct loop_it{
  size_t nd, i, iter;
  iterator it, beg;
  

  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
	typedef typename iterator_traits<iterator>::value_type value_type;
	typedef typename iterator_traits<iterator>::difference_type difference_type;
	typedef difference_type distance_type;	// retained
	typedef typename iterator_traits<iterator>::pointer pointer;
	typedef typename iterator_traits<iterator>::reference reference;


  loop_it(){}
  loop_it(iterator it_, size_t nd_, size_t iter_=0):it(it_), beg(it_), i(0), iter(iter_), nd(nd_){}
  
  /*dup_it(const dup_it &operator other):it(other.it),i(other.i), nd(other.nd){}
  dup_it& operator=(const dup_it &operator other){
    if(this!=&other)return *this=
  }*/

  value_type &operator*(){
    return *it;
  }
  loop_it &operator++(){
    i++;
    if(i<nd)++it;
    else{
      i=0;
      it=beg;
      iter++;
    }
    return *this;
  }
  loop_it operator++(int){ // postfix
    loop_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const loop_it &other) const {
    return (i==other.i && iter==other.iter);
  }
  bool operator!=(const loop_it &other) const {
    return !(*this==other);
  }

};


/// duplicator: duplicates each entry in the sequence nd times
template<class iterator, class checker_t=void_check_t>
struct dup_it{
  size_t nd, i;
  iterator it;

  typedef typename iterator_traits<iterator>::iterator_category iterator_category;
  typedef typename iterator_traits<iterator>::value_type value_type;
  typedef typename iterator_traits<iterator>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<iterator>::pointer pointer;
  typedef typename iterator_traits<iterator>::reference reference;

  dup_it(){}
  dup_it(iterator it_, size_t nd_, size_t i_=0):it(it_), i(i_), nd(nd_){}
  
  
  value_type &operator*(){
    return *it;
  }
  dup_it &operator++(){
    i++;
    if(i>=nd){
      ++it;
      i=0;
    }
    return *this;
  }
  dup_it operator++(int){ // postfix
    dup_it tmp=*this;
    ++*this;
    return tmp;
  }

  bool operator==(const dup_it &other) const {
    return (i==other.i && it==other.it);
  }
  bool operator!=(const dup_it &other) const {
    return !(*this==other);
  }

};


/// returns the first std::pair value of iterate pairs
template<class pair_it>
class first_it:public pair_it{
public:
  typedef typename iterator_traits<pair_it>::iterator_category iterator_category;
  typedef typename iterator_traits<pair_it>::value_type::first_type value_type;
  typedef typename iterator_traits<pair_it>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<pair_it>::pointer pointer;
  typedef typename iterator_traits<pair_it>::reference reference;

  first_it(){}
  
  first_it(const pair_it &base):pair_it(base){}
  
  value_type &operator*(){
    return ((pair_it *)this)->first;
  }  
};

/// returns the first pair value of iterate pairs
template<class pair_it>
class second_it:public pair_it{
public:
  typedef typename iterator_traits<pair_it>::iterator_category iterator_category;
  typedef typename iterator_traits<pair_it>::value_type::second_type value_type;
  typedef typename iterator_traits<pair_it>::difference_type difference_type;
  typedef difference_type distance_type;	// retained
  typedef typename iterator_traits<pair_it>::pointer pointer;
  typedef typename iterator_traits<pair_it>::reference reference;

  second_it(){}
  
  second_it(const pair_it &base):pair_it(base){}
  
  value_type &operator*(){
    return ((pair_it *)this)->second;
  }  
};


template <class T>
struct plusplus_t{
  T &operator()(T& v) const { v++; return v; }
};

template <class T>
struct minusminus_t{
  T &operator()(T& v) const { v--; return v; }
};

template <class T>
struct noop_t{
  T &operator()(T& v) const { return v; }
};

//e incremental iterator (to use together with index sets)
//e increment/decrement operators are applied to the value
//e returns the value itself when calling *it
template<class T, class incr_op=plusplus_t<T>, class decr_op=minusminus_t<T> >
class value_it{
protected:
  T v;
  incr_op incr;
  decr_op decr;
public:
  typedef random_access_iterator_tag iterator_category;
  typedef T value_type;
  typedef T difference_type;
  typedef difference_type distance_type;	// retained
  typedef T* pointer;
  typedef T& reference;

  value_it(){}
  value_it(T v_):v(v_){}
  
  T &operator*(){
    return v;
  }
  
  value_it &operator++(){
    //v++;
    incr(v);
    return *this;
  }

  value_it operator++(int){ // postfix
    value_it tmp=*this;
    ++*this;
    return tmp;
  }

  value_it &operator--(){
    //v--;
    decr(v);
    return *this;
  }

  // difference
  T operator-(const value_it & other) const {
    return v-other.v;
  }

  value_it operator--(int){ // postfix
    value_it tmp=*this;
    --*this;
    return tmp;
  }

  bool operator==(const value_it &other) const {
    return (v==other.v);
  }
  bool operator!=(const value_it &other) const {
    return !(*this==other);
  }

};

//e forward-iterating a random acces iterator via index set
//e increment, comparison operators  are applied to indexer
//e *it returns the value from rand_it by incrementing begin with the current index
template<class rand_it, class indexer_it>
class indexed_it{
  rand_it beg;
  indexer_it ind;
public:
  
public:
  typedef forward_iterator_tag iterator_category;
  typedef typename iterator_traits<rand_it>::value_type value_type;
  typedef typename iterator_traits<rand_it>::difference_type difference_type;
  typedef typename iterator_traits<rand_it>::distance_type distance_type;	// retained
  typedef typename iterator_traits<rand_it>::pointer pointer;
  typedef typename iterator_traits<rand_it>::reference reference;

  indexed_it(){}
  indexed_it(rand_it beg_, indexer_it ind_):beg(beg_),ind(ind_){}

  
  value_type &operator*(){
    return *(beg+*ind);
  }
  
  indexed_it &operator++(){
    ind++;
    return *this;
  }

  indexed_it operator++(int){ // postfix
    indexed_it tmp=*this;
    ++*this;
    return tmp;
  }

 
  bool operator==(const indexed_it &other) const {
    return (ind==other.ind);
  }
  bool operator!=(const indexed_it &other) const {
    return !(*this==other);
  }
 
};


/// Class for use in average_files
struct av_comm_t{
  string text;
  int flag;
};


/// aligned continuous memory allocation with allocation border of mult*sizeof(T) bytes
/// returns original pointer in orig_ptr which must be deleted by delete [] (T*)
template <class T>
inline T *aligned_alloc(size_t mult, size_t size, T **orig_ptr) {
  size_t nmult=size/mult;
  if(size%mult)
    nmult++;
  nmult++;
  *orig_ptr = new T[mult*nmult];
  
  unsigned long long pv=(unsigned long long)*orig_ptr, rest, div=mult*sizeof(T);
  rest=pv%div;
  if(rest)
    pv+=(div-rest);
  return reinterpret_cast<T*>(pv);

  //int mask = mult*sizeof(T)-1;
  
  //return (T *)(((size_t)((char *)*orig_ptr+mask))&(~mask));
}

/// avrerages the contents of space-separated text files containing numbers into the output
/// maximal number of columns: 100
/// maxinal line length: 10000 symbols
/// '#' symbols are treated as comments and transferred to output from the first valid file in sequence  
template < class filename_it>
int average_files(filename_it beg, filename_it end, const char *output){
  std::vector<av_comm_t> comments; 
  char buff[10000];
  double values[100];
  int ncol=0, nfiles=0; 
  std::vector< std::vector<double> > gtable;
 
  bool fix_comm=false;
  FILE *f=NULL;
  for(;beg!=end;++beg){
    std::vector< std::vector<double> > table(ncol);
    if(f)
      fclose(f);
    if(!fix_comm)
      comments.clear();
    f=fopen((*beg).c_str(),"rt");
    if(!f){
      LOGMSG(vblWARN,fmt("average_files: can't open the file '%s', skipped!",(*beg).c_str()),0);
      continue;
    }
    int line=0;
    size_t iline=0;
    bool skip=false;
    while(fgets(buff,10000,f)){
      ++line;
      char *comm=strstr(buff,"#");
      av_comm_t s;
      s.flag=-1;
      if(comm){
        s.text=comm;
        *comm=0;
      }
      if(!fix_comm)
        comments.push_back(s);
      
      int num=string_scan_t(buff,"%lf"," ",values,values+100);
      if(num<=0)// completely failed, next string
        continue;
      if(!ncol){
        ncol=num;
        table.resize(ncol);
      }
      else if(num<ncol){
        LOGMSG(vblWARN,fmt("average_files: unexpected number of entries at line %d in '%s', file skipped!",line,(*beg).c_str()),0);
        skip=true;
        break;
      }
      ++iline;
      for(int i=0;i<ncol;i++){
        if(table[i].size()<iline)
          table[i].resize(iline,values[i]);
        else 
          table[i][iline-1]+=values[i];
      }
      if(!fix_comm)
        comments.back().flag=1;
    }
    if(skip)
      continue;
    
    if(!gtable.size())
      gtable=table;
    else if(gtable[0].size()>table[0].size()){
      LOGMSG(vblWARN,fmt("average_files: number of lines in '%s' is too small, file skipped!",(*beg).c_str()),0);
      continue;
    }
    else{ // summing up
      for(size_t i=0;i<gtable.size();i++)
        for(size_t j=0;j<gtable[i].size();j++)
          gtable[i][j]+=table[i][j];
    }
    fix_comm=true; // fixing comments
    nfiles++;
  }
  if(f)
    fclose(f);
  if(!fix_comm)
    return LOGERR(0,fmt("average_files: no valid files with column data found!"),0);
  
  f=fopen(output,"wt");
  if(!f)
    return LOGERR(-1,fmt("average_files: can't open the file '%s' for writing",output),0);
  // retaining comment structure of the first valid file
  size_t iline=0;
  for(size_t i=0;i<comments.size();i++){
    if(comments[i].flag>0){
      for(size_t j=0;j<gtable.size();j++)
        fprintf(f,"%g ",gtable[j][iline]/nfiles);
      ++iline;
    }
    fprintf(f,"%s",comments[i].text.c_str());
    if(!comments[i].text.size())
      fprintf(f,"\n");
  }
  fclose(f);
  return 1;
}

///\en Writes a set of functions given by [\a fbeg, \a fend) into a multicolumn ASCII file
///    evaluated over a set of arguments [\a abeg, \a aend).
///    Arguments and functions must be compatible with the printf formats
///    supplied in \a arg_format and \a func_format respectively.
///    Arguments are printed in the first column, f(arg) values in the next columns.
///    Header is printed at the beginning,  tailer at end of file (with no extra \n after them).
///    Mode \a mode is supplied to fopen second argument (can be "wt" or "at").
///    \return >0 if OK.
template <class arg_it, class func_it>
int write_ascii(const char *fname, arg_it abeg, arg_it aend, 
               func_it fbeg, func_it fend, const char *header="",
               const char *arg_format="%g",const char *func_format="%g", const char *tab=" ",
               const char *mode="wt", const char *tailer=""){
  FILE *f=fopen(fname,mode);
  if(!f)
    return LOGERR(-1,fmt("write_ascii: can't open file '%s' intended for writing (mode: %s)!",fname,mode),0);
  fprintf(f,"%s",header);
  for(;abeg!=aend;++abeg){
    fprintf(f,arg_format,(*abeg));
    for(func_it fit=fbeg;fit!=fend;++fit){
      fprintf(f,tab);
      fprintf(f,func_format,(*fit)(*abeg));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"%s",tailer);
  fclose(f);
  return 1;
}


# if 0
//e incremental iterator reading values from the list
template<class T>
class va_it{
protected:
  T v;
  va_list argptr;
public:
  typedef forward_iterator_tag iterator_category;
  typedef T value_type;
  typedef int difference_type;
  typedef difference_type distance_type;	// retained
  typedef T* pointer;
  typedef T& reference;

  va_it(const T &first, ...){
    va_start(argptr, first);
    v = va_arg( argptr, T);
  }
  
  T operator*(){
    return v;
  }
  
  va_it &operator++(){
    v = va_arg( argptr, T);
    *this;
  }

  va_it operator++(int){ // postfix
    value_it tmp=*this;
    ++*this;
    return tmp;
  }

  //value_it &operator--(){
  //  //v--;
  //  decr(v);
  //  return *this;
 // }

  // difference
  //T operator-(const value_it & other) const {
  //  return v-other.v;
  //}

  va_it operator--(int){ // postfix
    value_it tmp=*this;
    --*this;
    return tmp;
  }

  bool operator==(const va_it &other) const {
    return (argptr==other.argptr);
  }
  bool operator!=(const va_it &other) const {
    return !(*this==other);
  }

};

# endif

# endif
