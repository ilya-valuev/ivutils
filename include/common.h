/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.24 $
 *   $Date: 2015/11/19 08:46:44 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/common.h,v 1.24 2015/11/19 08:46:44 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/common.h,v $
$Revision: 1.24 $
$Author: valuev $
$Date: 2015/11/19 08:46:44 $
*/
/*s****************************************************************************
 * $Log: common.h,v $
 * Revision 1.24  2015/11/19 08:46:44  valuev
 * common base for AWP
 *
 * Revision 1.23  2015/01/30 09:08:04  valuev
 * working on sweeps
 *
 * Revision 1.22  2012/01/17 06:00:07  valuev
 * neb workflow
 *
 * Revision 1.21  2009/08/22 12:45:35  morozov
 * Conditional compilation based on USE_MPI
 *
 * Revision 1.20  2009/08/21 18:25:03  morozov
 * Added basic MPI procedures
 *
 * Revision 1.19  2009/07/24 05:08:46  valuev
 * Sync with FDTD, added molecule setup
 *
 * Revision 1.18  2009/04/27 18:33:02  valuev
 * Fixed PBC in Coulomb potential: now supports mdPBC_EXTRA
 *
 * Revision 1.17  2009/03/04 09:52:45  valuev
 * corrected after sync with FDTD project
 *
 * Revision 1.16  2008/08/18 21:40:09  valuev
 * added Gurski-Krasko potential
 *
 * Revision 1.15  2008/07/23 16:42:08  valuev
 * Added AWPMD Monte-Carlo
 *
 * Revision 1.14  2008/06/03 15:23:37  valuev
 * Added rotation
 *
 * Revision 1.13  2008/04/22 22:26:55  valuev
 * working awpmd test for hydrogen molecules
 *
 * Revision 1.12  2008/04/22 12:44:17  valuev
 * made gcc 4.12 compilable
 *
 * Revision 1.11  2008/03/18 17:17:21  valuev
 * corrected PBC in microfield calculations
 *
 * Revision 1.10  2008/02/21 14:02:51  valuev
 * Added parametric methods
 *
 * Revision 1.9  2007/11/25 18:21:26  valuev
 * Added math switch to common
 *
 * Revision 1.8  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.11  2007/03/22 15:04:05  lesha
 * bfseek is added
 *
 * Revision 1.10  2007/03/20 17:58:09  lesha
 * accdiv moved to ifdef cplusplus
 *
 * Revision 1.9  2007/03/20 17:21:14  lesha
 * accdiv is added
 *
 * Revision 1.8  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.7  2006/12/25 11:53:15  valuev
 * Added dump for dispersive materials, fixed some errors
 *
 * Revision 1.6  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
 * Revision 1.5  2006/11/24 16:27:10  lesha
 * ifndef UNIX
 * define FU
 * endif
 *
 * is added
 *
 * Revision 1.4  2006/11/24 11:08:59  lesha
 * FU is #defined
 *
 * Revision 1.3  2006/10/27 20:41:01  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.5  2006/10/09 12:58:56  bogomolov
 * patch for compilling under unix
 *
 * Revision 1.4  2006/09/26 10:59:42  valuev
 * Added nonorthogonal TB (Menon-Subbaswamy)
 *
 * Revision 1.3  2006/03/14 10:32:17  valuev
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
# ifndef __COMMON_H
# define __COMMON_H

# ifndef _USE_MATH_DEFINES
# define _USE_MATH_DEFINES
# endif 
# include <math.h>

# include <stdio.h>
# include <stdlib.h>

/*# ifndef UNIX
# include <mem.h>
# else
# include <string.h>
# endif*/

# include <string.h>

# ifdef USE_MPI
# include <mpi.h>
# endif

#ifndef UNIX
# ifndef FU
# define FU
# endif
#endif

# ifndef NO_CMNMATH

/*
# ifndef min
# define min(x,y)  ((x)<(y)? (x): (y))
# endif

# ifndef max
# define max(x,y)  ((x)>(y)? (x): (y))
# endif
*/

# ifndef fmod
# define fmod(a,b)  ((a)-((long)((a)/(b))*(b)))
# endif

double log2(double x);

# endif

#define CMN_LOG2 0.693147180559945

typedef double realtype;
#define REALFRM "%lf" 



// some compilers don't define PI!
# ifndef M_PI
# define M_PI 3.1415926535897932385
# endif

# ifndef MAXINT
# define MAXINT 0x7fffffff;
# endif


# ifndef __BCC
# define random(num)  ((int)(((double)rand())/RAND_MAX*(num)))
# endif

//e returns normally distributed number with
//e variance delta and max. value of length
//e uses random1
double gaussrand(double length,double delta);

///en Returns 3 values in v as components of Maxwell-distributed
///   quantity with component-average norm sqauare of v_sq/2 (the average sum of 3 values is 3*v_sq/2).
void maxwellrand(double *v,double v_sq);


// returns maximum and minimum of realtype array
void minmax(realtype *arr, int n,int *imin, realtype *amin,
	    int *imax, realtype *amax);


# ifdef __cplusplus

# ifndef NO_BFSEEK

// for seek in files > 2 Gb
inline int bfseek(FILE *stream, long long offset, int whence) {
# ifdef UNIX
  return fseeko(stream, offset, whence);
# else
  return _fseeki64(stream, offset, whence);
# endif
}

# endif

# define VEC_ZERO3 1e-8

# if 0
// division accounting possible computational error
template <class num_t, class denum_t>
denum_t accdiv(const num_t x, const denum_t y) {

    denum_t id=x/y;
    denum_t fid=floor(id);
    if (id-fid<VEC_ZERO3)
      return fid;
    else if (id-fid>1.-VEC_ZERO3)
      return fid+1.;
    else
      return id;

/*
  int yexp;
  frexp(y, &yexp); // y=matissa*2^yexp; 0.5<=mantissa<1
  denum_t eps=numeric_limits<denum_t>::epsilon(); // possible error of 1
  denum_t dy=ldexp(eps, yexp-1); // possible error of dy
    
  denum_t fr=x/y;
  denum_t ifr=floor(fr+.5); // nearest integer
  denum_t diff=fr-ifr;
  denum_t err=fr/y*dy; // possible error of division
  
  return (fabs(diff)<2*err) ? ifr : fr;*/
}
# endif

/*
//e Verbosity level codes
enum VERB_LEVELS{
     VL_INFO3 =0x1,  //e< deep-level development info
     VL_INFO2 =0x2,
     VL_INFO1 =0x4,  //e< normal user info
     VL_INFO  =VL_INFO1|VL_INFO2,
     VL_MSG   =0x4,
     VL_WARN  =0x8,
     VL_ERR   =0x10,
     VL_NORMAL=VL_ERR|VL_WARN|VL_MSG,

     VL_ALL   =0xff,
     VL_ALARM =VL_ERR|VL_WARN,
     VL_NONE  =0
};*/

//e given a string with the opening 'bracket' before its beginning
//e position, returns the matching closing bracket index, or -1 if none found
//e if skip_cont is specified, it points to the bracket pairs around context to be skipped
//e for example if br_pair="()",  skip_cont="[]{}" and serach string is ("hjsdhk{[fhj(kk]}kk)"
int get_matching(const char *str,const char *br_pair, const char *skip_cont=NULL, const char *quotes=NULL, int from_back=0);


double vnorm2(double x1, double y1, double x2, double y2, double *nx=NULL, double *ny=NULL);
double vangle(double nx, double ny);


char *form(char *form,...);

extern "C" {

# endif

# ifdef __cplusplus
void serror(char *str,...);

extern void (*msg_error)(char *format, ...);
extern void (*fatal_error)(char *format, ...);
extern int (*eprintf)(const char *format,...);

extern int err_sts;
extern int err_status();


# else
void serror(char *str,...);



extern void (*msg_error)(char *format, ...);
extern void (*fatal_error)(char *format, ...);
extern int (*eprintf)(const char *format,...);

extern int err_sts;
extern int err_status();

# endif


typedef void (*err_func)(char *,...);
typedef int (*print_func)(const char *,...);
typedef int (*sask_func)(char *answers[],int *ind,char *format,...);

// asks the question using format and returns the
// answer choosen by user
int sask(char *answers[],int *ind,char *format,...);
int syes_no_cancel(char *format, ...);

extern sask_func ask_user;
extern print_func  ask_yes_no_cancel;

void shift_array(void *array,unsigned int top_arr,unsigned int start,int number);


/* exchanges size bytes between pointers a  and b  */
/* leaves size bytes allocated in memory ! */
void swapmem(void *a,void *b,size_t size);

extern realtype FZero_general;
# define Set_Zero_Level(a) {FZero_general=a;}

realtype lin_approx(realtype a,realtype f1,realtype b,realtype f2,realtype x);
long lin_approxl(long a,long f1,long b,long f2,long x);

FILE *Err_fopen(const char *file,const char *mode);
int fgetline(FILE *f,char *str,int len_tresh);

// skips end-of lines and
// lines beginning with comment
int fskip_comment(FILE *f,const char *comment);


# ifdef __cplusplus
}



// interpreting spec, example: 
// 10.22
// [10.,12.,0.1]
// returns 0 if bad format
// -num if error reading num's argument
// 1 if single par
// 3 if range
// sets corresponding values
int scan_range(char *str, realtype *a1, realtype *a2, realtype *step, char *delim=",");

// scans values from string into the buffer ptr
// the buffer is an array of the entries of size fieldsize
// the values are read according to format and separated by delimiter
// ends with the last value read or maxcount reached
// returns the number of fields read
int string_scan(const char *string,char *format,char *delim,void *ptr,size_t fieldsize,int maxcount);

typedef double (*dfuncp)(double);

dfuncp func_by_name(char *fname);


//e clears  array with optional integer index
void clear_arri(int n,double *vec, int *ind=0);

class NamedList{
 char **array;
 int nalloc;
public:
 int maxlen;
 int n;
 NamedList(int l=10){
  n=0;
  maxlen=l;
  nalloc=0;
 }
 ~NamedList();
 int search(char *name);
 int insert(char *name);
 int check(char *name){
  int i=search(name);
  if(i<0)return insert(name);
  return i;
 }
 char* operator[](int i){
  return array[i];
 }
};


class Pool{
 void **data;
 void *myptr;
public:

 int len;
 int limit;  // current allocated size
 int size1;  // first reserved amount
 int size_add;    // amount to be added to realloc data
 int ind; // current accessible size
 Pool(void **ptr=NULL,int fieldlen=sizeof(int)){
  if(!ptr)ptr=&myptr;
  data=ptr;
  *data=NULL;
  len=fieldlen;
  limit=0;
  ind=0;
  size1=10;
  size_add=10;
 }
 ~Pool(){
  if(limit && data==&myptr && *data)free(*data);
 }
 void clear(){
  if(limit && *data && data==&myptr)free(*data);
  limit=0;
  ind=0;
 }
 void reset(){
  ind=0;
 }   
 void allclear(){
  if(*data)free(*data);
  limit=0;
  ind=0;
 } 
 void set_adjust(int a){
  if(a)size_add=abs(size_add);
  else size_add=-abs(size_add);
 }
 int add(void *ptr,int n=1);
 void *operator[](int i){
   if(i>=ind)add(NULL,limit-i+1);
   return ((char*)*data)+i*len;
 }

 int set(int i, void *ptr){
   void *dtp=(*this)[i];
   memcpy(dtp,ptr,len);
   return ind;
 }

 void *get_ptr(){
   return *data;
 }
};


typedef Pool *PoolP;


// returns the coefficient d
// fav = f(-)*d+f(+)*(1-d)

realtype DxScale(realtype dx,realtype x,int &k);


enum FTYPES { XYDATA,YDATA};
enum BTYPES { OUT_ZERO,OUT_LAST};

typedef int (*TableMapFunc)(realtype *x, realtype *y);

class TableFunction{
private:
 void sort();
 static TableFunction* Fcur;
//public:
protected:
 char alloc;
 realtype *xx;
 realtype *yy;
public:
 enum FTYPES ftype;
  // public:
 friend realtype TabFunc(realtype x);
 friend void SetFunc(TableFunction *f);

 int n;
 enum BTYPES btype;
 realtype xstart,xend;
 TableFunction(){ xx=yy=NULL; alloc=0;n=0; }
 TableFunction(realtype x1, realtype x2, int dim);
 TableFunction(int dim,realtype *arr);
 TableFunction(int dim,realtype *arrx,realtype *arry);
 TableFunction(char *file,int col1,int col2);
 TableFunction(FILE *f, int col1, int col2);
 void write(char *file,realtype yscale=1.,TableMapFunc mfunc=NULL);
 virtual ~TableFunction(){
  if((alloc&0x01))delete [] yy;
  if((alloc&0x02))delete [] xx;
  xx=yy=NULL;
  n=0;
 }

 // trapecial functions
 virtual realtype integral();
 // returns the first x where integral reaches val
 // if not, returns xend
 // works ONLY for positively defined functions!!!
 virtual realtype integral_reaches(realtype val);

 // integrates current function, result is intergral(x)
 // scales the value with the given scale factor
 virtual TableFunction& integrate(realtype scale=1.);

 /* avoid using this function */
 virtual TableFunction& operator=(const TableFunction& F){
   if((alloc&0x01))delete[] yy;
   if((alloc&0x02))delete[] xx;
   //eprintf("ehe\n");
   if(F.alloc){
     eprintf("TableFunction: equating allocated item!\n");
   }
   memcpy((void*)this,(void*)&F,sizeof(TableFunction));
   alloc=0;
   return *this;
 }

 TableFunction(TableFunction &F){
   //eprintf("oo\n");
   if(F.alloc){
     eprintf("TableFunction: copying allocated item!\n");
   }
   memcpy((void*)this,(void*)&F,sizeof(TableFunction));
   alloc=0;
   F.alloc=0;
 }


  // copy creation
 virtual TableFunction& operator<<(const TableFunction &F);



 virtual TableFunction& operator*=(TableFunction &F);
 virtual TableFunction& operator*=(realtype f(realtype));
 virtual TableFunction& operator*=(realtype c);
 virtual TableFunction& operator/=(TableFunction &F);
 virtual TableFunction& operator+=(TableFunction &F);
 virtual TableFunction& operator+=(realtype f(realtype));
 virtual TableFunction& operator+=(realtype c);
 virtual TableFunction& operator-=(TableFunction &F);
 virtual TableFunction& operator-=(realtype f(realtype));
 virtual TableFunction& operator-=(realtype c);

 virtual TableFunction& operation(TableFunction &F,realtype fop(realtype,realtype));  

 realtype x(int i){
  if(ftype==YDATA)return xstart+(xend-xstart)*i/(realtype)(n-1);
  return xx[i];
 }
 realtype &y(int i){ return yy[i];}
 realtype xscale(realtype x1,realtype x2);
 int   nsteps,stpk,stpi,insert_stat;
 realtype  stpx,yold;
 void  insert_begin(int steps);
 realtype insert_next(realtype y);
 
 virtual realtype operator()(realtype x);
 virtual realtype der(realtype x);

 virtual realtype Dx(realtype x);
 virtual realtype Dy(realtype x); 
};

extern  realtype TabFunc(realtype x);
extern  void SetFunc(TableFunction *f);


// reads one- or two- column file into the arrays
// allocates memory in yy and (if needed) xx.
// Returns the number of entries
int ReadTable(char *file,realtype* &xx,realtype* &yy,int col1,int col2);

int ReadTableDirect(FILE *f,realtype* &xx,realtype* &yy,int col1,int col2,int ws);


// Sets array in order according to 'order' list
void SetInOrder(int *order,int n, void *array, size_t size);

// compare function to sort realtype array indexed by integer
// order list by qsort
// if CmpAction==1 sorts in ascending order
// if CmpAction==-1 sorts in descending order
extern realtype *CmpArray;
extern int CmpAction;
int CompIndexed(const void *c1,const void *c2);


class SplineFunction:public TableFunction{
 realtype *aa;
 realtype *bb;
 realtype *cc;
 void makespl();
 void allocspl();
public:
 SplineFunction(int dim,realtype *arr);
 SplineFunction(int dim,realtype *arrx,realtype *arry);
 SplineFunction(char *file);
 ~SplineFunction(){
  if((alloc&0x04)){
   delete [] aa;
   delete [] bb;
   delete [] cc;
  }
 }
 realtype operator()(realtype x);
 realtype der(realtype x);
 realtype xscale(realtype x1,realtype x2);
};

# endif

// truncates spaces (' ') from the end and beginning of a string
// modifies its argument!
char *trunc_spaces(char *);
// the same based on isspace function
char *trunc_spaces2(char *);

// escapes characters like tabs and CR in a string
// modifies its argument!
char *escape_chars(char *str);
char *deescape_chars(char *str);

char *set_extension(char *result,char *fname,char *ext);
char *get_extension(char *fname);
char *get_filename(char *fname);
char *get_directory(char *fname);
char *get_basename(char *fname);
int *GetList(int *n,char *str);

# ifdef __cplusplus

class cList{
 int *list;
 int pc;
 int ind;
 int vald;
 void copy_data(const cList &l);
public:
 static int stop_on_error;
 int n;
 cList(const char *str){
  char *buf=new char[strlen(str)+1];
  strcpy(buf,str);
  //printf("Creating (str) ...\n");
  list=GetList(&n,buf);
  delete [] buf;
  if(!list){
   if(stop_on_error)msg_error("Can't allocate list: %s\n",str);
   n=0;
  }
  pc=ind=0;
  if(n>0)vald=1;
  else vald=0;
 }
 cList(){
   list=NULL;
   pc=ind=vald=0;
   n=0;
   //printf("Creating (void)...\n");
 }
 cList(int start, int count){
   //printf("Creating (range)...\n");
  list=(int *)malloc(3*sizeof(int));
  if(!list){
   if(stop_on_error)msg_error("cList: error allocating memory.\n");
   n=0;
   pc=ind=vald=0;
   return;
  }
  list[0]=count;
  list[1]=start;
  list[2]=0;
  pc=ind=0;
  n=count;
  if(count>0)vald=1;
  else vald=0;
 }
 int current(){
  if(vald)return list[2*pc+1]+ind;
  else return -1;
 }
 int next();
 int step(){
  int cc=current();
  next();
  return cc;
 }
 int valid(){
  return vald;
 }
 int rewind(){
  pc=ind=0;
  if(n>0)vald=1;
  return 0;
 }
 // adjusts current entry to range and moves to the next entry
 // returns :
 // 0- no adjustment needed, 0x1|0x2 -- left| right boundary adjusted
 // 0x4 -- entry does not fit and disabled, -1 -- end of list
 int adjust_next(int i1, int i2);
 // adjusts the whole list to range
 int adjust(int i1, int i2);
 int current_range(int &i1, int &i2);

 // returns the number of occurences of index is in the list
 int is_in(int i);
 cList &operator=(const cList &);
 cList(const cList &l);
   
 ~cList();
 /*if(list)free(list);
 //}*/
};

class SymmMatr{
  int managed;
  double *arr;
public:
  
  unsigned size;

  SymmMatr():managed(0){
    size=0;
    arr=NULL;
  }

  int init(unsigned n){
    size=n;
    if(!managed)return 1;
    if(arr) delete [] arr;
    arr=new double[(long)n*(n+1)/2];
    if(!arr){
      msg_error("SymmMatr: memory allocation error.\n");
      return 0;
    }
    return 1;
  } 

  SymmMatr(unsigned n):arr(NULL),managed(1){ 
    init(n);
  }

  SymmMatr(unsigned n, double *ptr):arr(ptr),managed(0){ 
    init(n);
  }

  void Set(double val);

  void SetDiag(double val);

  virtual long GetDataPtr(double **ptr){
    *ptr=arr;
    return (long)size*(size+1)/2;
  }

  ~SymmMatr(){ 
    if(managed && arr)delete[] arr;
  }

  double &operator()(int i,int j){
    if(i>=j)return arr[(long)(2*size-j-1)*(long)j/2+i];
    else return arr[(long)(2*size-i-1)*(long)i/2+j];
  }
};

// for compatibility with the wrong syntax
typedef SymmMatr SimmMatr;


# endif


void Exit_wait(int val);


# ifndef ZERO
 # define ZERO ISZERO
# endif

# define ISZERO(a)           ( fabs(a)-FZero_general < 0 ? 1 : 0 )
# define SIGN(a)           ( fabs(a)-FZero_general < 0 ? 0 : ( a< 0 ? -1 : 1) )
# define fmin(a,b) ((a)>(b) ? (b):(a))
# define fmax(a,b) ((a)>(b) ? (a):(b))
# define fsign(a) ((a)>0 ? 1:-1)
# define Err_malloc(a,b) { (a)=malloc(b);if(!(a))serror("Memory  allocation error.");}

# define Err_farmalloc(a,b) { (a)=farmalloc(b);if(!(a))serror("Far Memory  allocation error.");}

# define swap_var(a,b,type) { type t;t=(a);(a)=(b),(b)=(t); }

//# pragma argsused
int no_output(const char *format, ...);

// returns the maximal common denominator
int com_div(int a, int b);

# ifdef __cplusplus 
long com_divl(long a, long b);

int cmp_realtypec(const void *c1,const void *c2);


# endif

long  com_divl(int n, long *a);



# define IN_FRAME(x,xmin,xmax) ((x)>=(xmin) && (x) <= (xmax))




//e returns random [-0.5,0.5]
extern int rand_init;
double random1(void);

# ifdef UNIX // define random generator

// int random_r(int);

//# define random random_r
# endif

# ifdef FU

#include <stdarg.h>

# ifdef __cplusplus
extern "C" {
# endif 

int vsscanf(const char *,const char *, va_list);

# ifdef __cplusplus
}
# endif


# endif

#endif



