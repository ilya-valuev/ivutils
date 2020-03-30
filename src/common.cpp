/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.20 $
 *   $Date: 2015/11/19 08:46:44 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/common.cpp,v 1.20 2015/11/19 08:46:44 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/common.cpp,v $
$Revision: 1.20 $
$Author: valuev $
$Date: 2015/11/19 08:46:44 $
*/
/*s****************************************************************************
 * $Log: common.cpp,v $
 * Revision 1.20  2015/11/19 08:46:44  valuev
 * common base for AWP
 *
 * Revision 1.19  2009/08/22 21:28:34  morozov
 * Corrected MPI version
 *
 * Revision 1.18  2009/08/22 12:45:35  morozov
 * Conditional compilation based on USE_MPI
 *
 * Revision 1.17  2009/08/21 18:25:03  morozov
 * Added basic MPI procedures
 *
 * Revision 1.16  2008/09/22 20:43:51  valuev
 * added new repulsive C+O, fixed Gurski potential
 *
 * Revision 1.15  2008/08/18 21:40:09  valuev
 * added Gurski-Krasko potential
 *
 * Revision 1.14  2008/04/22 22:26:55  valuev
 * working awpmd test for hydrogen molecules
 *
 * Revision 1.12  2008/04/21 23:13:44  valuev
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
 * Revision 1.9  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.8  2007/01/22 21:27:50  valuev
 * Added domain decomposition, fixed DetectorSet (transform)
 *
 * Revision 1.7  2006/11/27 09:39:13  valuev
 * Added separators and plane mode to detectors
 *
 * Revision 1.6  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
 * Revision 1.5  2006/11/24 11:07:18  lesha
 * Number of possible elements in vsscanf is increased from 10 to 15
 *
 * Revision 1.4  2006/10/27 20:41:02  valuev
 * Added detectors sceleton. Updated some of ivutils from MD project.
 *
 * Revision 1.3  2006/10/19 21:42:18  lesha
 * log is removed
 * 6 -> 12 arguments is corrected
 *
 * Revision 1.2  2006/08/08 13:02:04  valuev
 * Added geometry
 *
 * Revision 1.3  2006/03/14 10:32:17  valuev
 * Added SetControls support for many components,
 * improved tcpenfine, added GRASP interface
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:16  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/

# ifndef __COMMON__
# define __COMMON__  

# ifdef USE_STDAFX
# include "stdafx.h"
# endif

# ifdef UNIX
# define cdecl
# endif

# include <stdio.h>
# include <stdarg.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <limits>
//# include <alloc.h>
# include "common.h"




# ifndef NO_CMNMATH

double log2(double x){
  return log(x)/CMN_LOG2;
}

# endif


double vnorm2(double x1, double y1, double x2, double y2, double *nx, double *ny){
  double x=x2-x1;
  double y=y2-y1;
  double val=sqrt(x*x+y*y);
  if(val<1e-32)return 0.;
  if(nx)*nx=x/val;
  if(ny)*ny=y/val;
  return val;
}

double vangle(double nx, double ny){
  double theta;
  if(fabs(nx)<1e-10){
    if(ny>0)theta=M_PI/2;
    else theta=-M_PI/2;
  }
  else{
    theta=atan(ny/nx);
    if(nx<0)theta+=M_PI;
  }
  return theta;
}

double gaussrand(double length,double delta){
  double tmp,val,x;
  do{
    x=2*length*random1();
    val=exp(-x*x/delta);
    tmp=0.5+random1();
    //printf("" REALFRM "  " REALFRM "\n",tmp,val);
  }while(tmp>val);
  return x;
}

void maxwellrand(double *v,double v_sq){
  const double epsilon = std::numeric_limits<double>::min();
  double xi[4];
  int i;
  for(i=0;i<4;i++){
    do{
      xi[i]=  ((double)rand())/RAND_MAX;  //0.5+random1();
    }while(i<2 && xi[i]<=epsilon); // otherwise log will be infinite
  }

  double a1=sqrt(fabs(-2.*v_sq*log(fabs(xi[0]))));
  /*if(!_finite(a1)){
    a1 = std::numeric_limits<double>::max()/2.;
    printf("%g %g\n", xi[0], epsilon);
  }*/
  double a2=sqrt(fabs(-2.*v_sq*log(fabs(xi[1]))));
  /*if(!_finite(a2)){
    a2 = std::numeric_limits<double>::max()/2.;
    printf("%g %g\n", xi[1], epsilon);
  }*/

  v[0]=a1*cos(2*M_PI*xi[2]);
  v[1]=a1*sin(2*M_PI*xi[2]);
  v[2]=a2*cos(2*M_PI*xi[3]);
  
}


void clear_arri(int n,double *vec, int *ind){
  if(!vec)return;
  int i, c;
  for(i=0;i<n;i++){
    if(ind)c=ind[i];
    else c=i;
    vec[c]=0.;
  }
}



struct {
char *name;
dfuncp func;
} DFuncTable[]={
  { "sin", sin },
  { "cos", cos },
  { "sqrt", sqrt },
  { "fabs", fabs },
  { "exp", exp },
  { 0,0}
};


dfuncp func_by_name(char *fname){
  int i=0;
  
  while(DFuncTable[i].name){
    if(!strcmp(fname,DFuncTable[i].name)){
      break;
    }
    i++;
  }
  return DFuncTable[i].func;
}





void minmax(realtype *arr, int n,int *imin, realtype *amin,
	    int *imax, realtype *amax){
  int i;
  if(n<=0)return;
  *imin=*imax=0;
  *amin=*amax=arr[0];
  for(i=1;i<n;i++){
    if(arr[i]<*amin){
      *amin=arr[i];
      *imin=i;
    }
    if(arr[i]>*amax){
      *amax=arr[i];
      *imax=i;
    }
  }
}


int string_scan(const char *str,char *format,char *delim,void *ptr, size_t fieldsize, int maxcount){
  char *buf=new char[strlen(str)+1];
  strcpy(buf,str);
  char *cptr=(char *)ptr;
  char *c1=strtok(buf,delim);
  int c=0;
  while(c1 && c<maxcount){
    int res=sscanf(c1,format,cptr);
    if(res>0){
      c++;
      cptr+=fieldsize;
    }
    c1=strtok(NULL,delim);
  }
  delete [] buf;
  return c;
}


int scan_range(char *range, realtype *a1, realtype *a2, realtype *step, char const *delim){
  char *c1, *c2;
  char str[500];
 
  strncpy(str,range,500);
 
  c1=strstr(str,"[");
  if(c1){ // assuming range
    c2=strstr(str,"]");
    if(!c2)return 0;// invalid spec
    *c2=0;
    c1=strtok(c1+1,delim);
    if(c1){
      if(!sscanf(c1,"" REALFRM "",a1))return -1;
    }
    else return 0;
    c1=strtok(NULL,delim);
    if(c1){
      if(!sscanf(c1,"" REALFRM "",a2))return -2;
    }
    else return 0;
    c1=strtok(NULL,delim);
    if(c1){
      if(!sscanf(c1,"" REALFRM "",step))return -3;
    }
    else return 0;
    return 3;
  }
  else{
    if(!sscanf(str,"" REALFRM "",a1))return -1;
  }
  return 1;
}



int err_sts=0;

int sask(char *answers[],int *ind,char *format,...){
  char string[256], antw[256];
  va_list vl;
  va_start(vl,format);
  vsprintf(string,format,vl);
  strcat(string," (");
  int i=0;
  while(answers[i]){
    strcat(string,answers[i]);
    strcat(string,"/");
    i++;
  }
  string[strlen(string)-1]=')';

  do{
    printf("%s\n",string);
    scanf("%s",antw);
    i=0;
    while(answers[i]){
      if(strstr(answers[i],antw))break;
      i++;
    }
  }while(!answers[i]);
  return ind[i];
}

int syes_no_cancel(char *format,...){
  char string[256];
  va_list vl;
  va_start(vl,format);
  vsprintf(string,format,vl);

  char *answers[]={"y","n","cancel",NULL};
  int ind[]={1,0,-1};
  return sask(answers,ind,string);
}


# ifdef __cplusplus

char *form(char *form,...){
 static char string[256];
 va_list vl;
 va_start(vl,form);
 vsprintf(string,form,vl);
 return string;
}

int press_key=0;
void Exit_wait(int val){
  press_key=val;
}

void serror(const char *format,...){
 va_list vl;
 va_start(vl,format);
 vprintf(format,vl);
 if(press_key){
  printf("Press ENTER to exit...\n");
  getchar();
 }
 exit(13);
}


int err_status(){
  return err_sts;
}


// generic dialog interface
void (*msg_error)(const char *format, ...)= serror;
void (*fatal_error)(const char *format, ...)= serror;
# ifdef USE_MPI
int master_printf(const char *format,...) {
  va_list vl;
  va_start(vl,format);

  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  if(irank != 0) return 0;
  return vprintf(format,vl);
}
int (*eprintf)(const char *format,...)= master_printf;
# else
int (*eprintf)(const char *format,...)= printf;
# endif
sask_func   ask_user   = sask;
print_func  ask_yes_no_cancel = (print_func)syes_no_cancel;


/*
void serror(int errno,char *comments=" "){
 switch(errno){
  case MEM:
   serror("%sMemory allocation error.\n",comments);
  case FILE:
   serror("Can't open file %s*/
# else

void cdecl serror(char *str){
 printf("%s\n",str);
 exit(0);
}
# endif

void shift_array(void *array,unsigned int top_arr,unsigned int start,int number){
 char *arr=(char *)array;
 unsigned int i;
 if(number>0){
  for(i=top_arr-1;i<=start;i--){
   arr[i+number]=arr[i];
  }
 }
 else{
  number=-number;
  for(i=start;i<top_arr-number;i++){
   arr[i]=arr[i+number];
  }
 }
}

realtype FZero_general=1e-18f;

FILE* Err_fopen(const char *file,const char *mode){
 char str[30];
 FILE *f;
 f=fopen(file,mode);
 if(!f){
  sprintf(str,"Can't open file %s.\n",file);
  fatal_error(str);
 }
 return f;
}

realtype cdecl lin_approx(realtype a,realtype f1,realtype b,realtype f2,realtype x){
 if(a==b)return f1;
 return f1+(f2-f1)*(x-a)/(b-a);
}

long cdecl lin_approxl(long a,long f1,long b,long f2,long x){
 f1=f1<<16;
 f2=f2<<16;
 if(a==b)return f1;
 return f1+(f2-f1)*(x-a)/(b-a);
}


/* exchanges size bytes between pointers a  and b  */

void cdecl swapmem(void *a,void *b,size_t size){
 static void *c;
 unsigned i;

 Err_malloc(c,size);
 for(i=0;i<size;i++)((char *)c)[i]=((char *)a)[i];
 for(i=0;i<size;i++)((char *)a)[i]=((char *)b)[i];
 for(i=0;i<size;i++)((char *)b)[i]=((char *)c)[i];
 free(c);
}

int fgetline(FILE *f,char *str,int len_tresh=0){
 char ch;
 int i=0;
 ch=(char)fgetc(f);
 while(!feof(f)){
  if(ch=='\n')break;
  str[i++]=ch;
  if(i==len_tresh)break;
  ch=(char)fgetc(f);
 }
 str[i]=0;
 return i;
}

// skips end-of lines and
// lines beginning with comment
int fskip_comment(FILE *f,const char *comment){
  int c, i=0;
  while(!feof(f)){
    do{
      c=fgetc(f);
    }while(c=='\n' && !feof(f));
    if(c==comment[0]){ // first char found
      i=1;
      while(comment[i] && !feof(f)){ // checking others
        c=fgetc(f);
        if(comment[i]!=c){ // failed, must ungetc
          ungetc(c,f);
          do{
            i--;
            ungetc(comment[i],f);
          }while(i);
          return 0;
        }
        i++;
      }
      // others
      // skipping the next line
      while(!feof(f)){
        c=fgetc(f);
        if(c=='\n')break;
      };
      continue;
    }
    else ungetc(c,f);
    return 0;
  };
  return 0;
}

char *trunc_spaces(char *str){
 int i=0,j;
 if(!str)return str;
 while(str[i]==' ' /*|| str[i]=='\8'*/)i++;
 j=(int)strlen(str)-1;
 while(j>=0 && (str[j]==' '/* || str[j]=='\8'*/))j--;
 str[j+1]=0;
 return str+i;
}

char *trunc_spaces2(char *str){
 int i=0,j;
 if(!str)return str;
 while(isspace(str[i]))i++;
 j=(int)strlen(str)-1;
 while(j>=0 && (isspace(str[j])) )j--;
 str[j+1]=0;
 return str+i;
}

char _escape_tbl[]="t\tn\n8:";


char *escape_chars(char *str){
 char *buff=new char[2*strlen(str)+1];
 char *ii=str, *io=buff, *j;
 char t;
 while(*ii){
   t=*ii;
   j=_escape_tbl+1;
   while(*j){
     if(*j==t){
       *io++ ='\\';
       t=*(j-1);
       break;
     }
     j+=2;
   }
   *io++=t;
   ii++;
 }
 *io++=0;
 strcpy(str,buff);
 delete [] buff;
 return str;
}

char *deescape_chars(char *str){
 char *buff=new char[strlen(str)+1];
 int ii=0,io=0,j;
 char c,t,s;
 while(str[ii]){
   t=str[ii];
   if(t=='\\'){
     s=str[ii+1];
     j=0;
     while(_escape_tbl[j]){
       c=_escape_tbl[j];
       if(c==s){
         t=_escape_tbl[j+1];
         ii++;
         break;
       }
       j+=2;
     }
   }
   buff[io++]=t;
   ii++;
 }
 buff[io++]=0;
 strcpy(str,buff);
 delete [] buff;
 return str;
}




char *set_extension(char *result,char *fname,char *ext){
 char arr[200],*s;
 strcpy(arr,fname);
 s=strtok(arr,".");
 int l=strlen(s);
 s[l]='.';
 strcpy(s+l+1,ext);      
 strcpy(result,s);
 return result;
}

char *get_extension(char *fname){
 char *ext=strstr(fname, ".");
 if(ext==NULL)return NULL;
 else return ext+1;
}

char *get_filename(char *fname){
 static char result[250];
 int i=strlen(fname);
 for(;i>=0;i--){
  if(fname[i]=='\\' || fname[i]=='/')break;
 }
 if(i!=0)i++;
 strncpy(result,fname+i,250);
 return result;
}


char *get_basename(char *fname){
 static char result[250];
 int i=strlen(fname);
 for(;i>=0;i--){
  if(fname[i]=='\\' || fname[i]=='/')break;
 }
 if(i!=0)i++;
 
 char *ext=strstr(fname, ".");
 if(ext)*ext=0;
 strncpy(result,fname+i,250);
 if(ext)*ext='.';
 return result;
}





char *get_directory(char *fname){
 static char result[250];
 int i=strlen(fname);
 for(;i>=0;i--){
  if(fname[i]=='\\' || fname[i]=='/')break;
 }
 if(i!=0){
  i++;
  strncpy(result,fname,i);
 }
 else{
  result[0]=0;
 }
 return result;
}






NamedList::~NamedList(){
 if(n!=0){
  for(int i=0;i<n;i++){
   delete array[i];
  }
  free(array);
 }
}

int NamedList::search(char *name){
 if(n<=0)return -1;
 for(int i=0;i<n;i++){
  if(strcmp(array[i],name)==0)return i;
 }
 return -1;   
}

int NamedList::insert(char *name){
 if(n==0){
  nalloc=10;
  array=(char **)malloc(nalloc*sizeof(char *));
  if(!array)fatal_error("Error allocating memory for NamedList.\n");
 }
 if(n>=nalloc){
  nalloc+=maxlen;
  array=(char **)realloc(array,(nalloc)*sizeof(char *));
  if(!array)fatal_error("Error reallocating memory for NamedList.\n");
 }
 array[n]= new char[strlen(name)+1];
 if(!array[n])fatal_error("Memory allocation error in NamedList.\n");
 strcpy(array[n++],name);
 return n-1;
}



int Pool::add(void *ptr,int n){
 while(ind+n>=limit){
  if(limit==0){
   limit=size1;
   *data=malloc(limit*len);
   if(!(*data))fatal_error("Pool: allocation error.\n");
  }
  else{
   if(size_add<=0)fatal_error("Pool: too many entries (%d)\n",limit+1);
   limit+=size_add;
   *data=realloc(*data,limit*len);
   if(!(*data))fatal_error("Pool: reallocation error.\n");
  }
 }
 if(ptr)memcpy((char *)*data+ind*len,ptr,n*len);
 ind+=n;
 return ind;        
}



int cList::stop_on_error=0;

cList::cList(const cList &l){
  copy_data(l);
}

void cList::copy_data(const cList &l){
 ind=l.ind;
 vald=l.vald;
 pc=l.pc;
 n=l.n;
 if(!l.list)return;
 int i=0;
 while(l.list[2*i]!=0)i++;
 list=(int *)malloc((2*i+1)*sizeof(int));
 if(!list)fatal_error("Memory allocation error copying list..\n");

 i=0;
 while(l.list[2*i]!=0){
  list[2*i]=l.list[2*i];
  list[2*i+1]=l.list[2*i+1];
  i++;
 }
 list[2*i]=0;
}

// returns the list of positive numbers from a string and the number
// of entries in *num
// format: n1,first1,n2,first2,...,0
cList &cList::operator=(const cList &l){
  if(this!=&l){
    if(list)free(list);
    list=NULL;
    copy_data(l);
  }
  return *this;
}

 


cList::~cList(){
  //printf("Deleting...\n");
 if(list)free(list);
 ind=0;
}

int cList::next(){
 if(!vald)return -1;
 ind++;
 while(ind>=list[2*pc]){
  pc++;
  ind=0;
  if(list[2*pc]==0){
   vald=0;
   pc=0;
   return -1;
  }
 }
 return list[2*pc+1]+ind;
}


int cList::adjust_next(int i1, int i2){
  int retval=0;
  if(!vald)return -1;

  if(list[2*pc+1]<i1){
    n-=list[2*pc];
    list[2*pc]-=i1-list[2*pc+1];
    list[2*pc+1]=i1;
    retval|=0x1;
  }
  if(list[2*pc+1]>i2){
    n-=list[2*pc];
    list[2*pc]=-1;
    retval=0x4;
  }
  else if(list[2*pc]+list[2*pc+1]-1>i2){
    n-=list[2*pc];
    list[2*pc]=i2-list[2*pc+1]+1;
    retval|=0x2;
  }
  if(retval&0x3){
   if(list[2*pc]>0)n+=list[2*pc];
   else{
    list[2*pc]=-1;
    retval=0x4;
   }
  }
  pc++;
  while(list[2*pc]<0)pc++;
  if(list[2*pc]==0)vald=0;
 
  return retval;
}


int cList::adjust(int i1, int i2){
  rewind();
  while(adjust_next(i1,i2)>=0);
  rewind();
  return 1;
}

int cList::current_range(int &i1, int &i2){
  if(vald){
   i1=list[2*pc+1];
   i2=i1+list[2*pc]-1;
  }
  return vald;
}

int cList::is_in(int index){
 if(n<=0)return 0;
 int i=0;
 int c=0;
 while(list[2*i]!=0){
  if(index>=list[2*i+1] && index<list[2*i+1]+list[2*i])c++;
  i++;
 }
 return c;
}




int *GetList(int *num,char *str){
 int *list=NULL,n,n1,n2;
 char *str1,*p;
 Pool pl((void**)&list);
 str1=strtok(str,",");
 *num=0;
 if(!str1)return list;

 do{
  p=strstr(str1,"-");
  if(p){
   *p=0;
   if(sscanf(str1,"%d",&n1)!=1 || sscanf(p+1,"%d",&n2)!=1){
    if(list)free(list);
    return NULL;
   }
   n=abs(n2-n1)+1;
   pl.add(&n);
   if(n1>n2)pl.add(&n2);
   else pl.add(&n1);

   *p='-';
   (*num)+=n;
  }
  else{
   if(sscanf(str1,"%d",&n)!=1){
    if(list)free(list);
    return NULL;
   }
   n1=1;
   pl.add(&n1);
   pl.add(&n);
   (*num)++;
  }
  str1=strtok(NULL,",");
 }while(str1);
 n=0;
 pl.add(&n);
 return list;
}

// reads the table from a file
// if whole_scan==1 reads only till the first failure
// col1 -- x column number
// col2 -- y column number
// if col1<=0, reads only y, x is assumed to be [0,1.] range scale
// if col2<=0 --> autodetect, first two file columns are used

int ReadTableDirect(FILE *InFile,realtype* &xx,realtype* &yy,int col1, int col2, int whole_scan=1){
 
 char str[1000];
 realtype x,y;                                           
 int Count=0,InType=2,cur1,cur2,i,read_failed;
 char frm1[250]="", frm2[250]="";
 
 if(col2<=0){
   InType=0;
   col1=1;
   col2=2;
 }
 else if(col1<=0){
   col1=col2;
   col2=-1;
   InType=1;
 }
 
   
 for(i=0;i<col1-1;i++)strcat(frm1,"%*f ");
 strcat(frm1,"" REALFRM "");
 for(i=0;i<col2-1;i++)strcat(frm2,"%*f ");
 if(col2>=0)strcat(frm2,"" REALFRM "");

 //printf("frm: (%s, %s)\n",frm1,frm2);

 long initpos=ftell(InFile);
 i=0;
 do{
  fgetline(InFile,str,1000);

  cur1=sscanf(str,frm1,&x);
  if(cur1==EOF)cur1=0;
  cur2=sscanf(str,frm2,&y);
  if(cur2==EOF)cur2=0;
  //printf("C: %s <- %d %d\n",str,cur1,cur2);

  if(cur1 && InType<=(cur1+cur2)){
   if(!InType)InType=cur1+cur2;

   if(i){ //writing
    if(InType==1)yy[Count]=x;
    else{      
      xx[Count]=x;
      yy[Count]=y;
     
    }
   }
   Count++; //counting
   read_failed=0;
  }
  else read_failed=1;
  if(feof(InFile) || ((!whole_scan) && read_failed)){
   if(Count==0 || i!=0)break;

   if(yy) delete [] yy;
   yy=new realtype[Count];
   if(!yy)fatal_error("ReadTable: memory allocation error.\n");
   if(InType==2){
     if(xx) delete [] xx;
     xx=new realtype[Count];
     if(!xx)fatal_error("ReadTable: memory allocation error.\n");
   }
   i++;
   fseek(InFile,initpos,SEEK_SET);
   Count=0;
  }
 }while(i<2);

 return Count;
}

int ReadTable(char *file,realtype* &xx,realtype* &yy, int col1=-1, int col2=-1){
 FILE *InFile;
 InFile=Err_fopen(file,"rt");
 int Count=ReadTableDirect(InFile,xx,yy,col1,col2,1);
 fclose(InFile);
 return Count;
}


// returns the coefficient d
// fav = f(-)*d+f(+)*(1-d)

realtype DxScale(realtype dx,realtype x,int &k){
 k=(int)(x/dx);
 realtype r=x-dx*k;
 if(r<0)r=-r, k--;
 return r/dx;
}

TableFunction::TableFunction(int dim,realtype *arr=NULL){
 xx=NULL;
 yy=arr;
 n=dim;
 if(n<=0)fatal_error("TableFunction: Invalid array dimension: %d\n",n);
 alloc=0;
 if(yy==NULL){
   yy=new realtype[dim];
   if(!yy)fatal_error("TableFunction: memory allocation error.\n");
   alloc=0x1;
   for(int i=0;i<dim;i++)yy[i]=0.;
 }
 xstart=0.;
 xend=1.;
 ftype=YDATA;
 btype=OUT_ZERO;
}

TableFunction::TableFunction(realtype x1, realtype x2,int dim){
  TableFunction a(dim, NULL);
  *this=a;
  a.alloc=0;
  xstart=x1;
  xend=x2;
}

/*
A kak zhe slovari? Neuzheli
ty povezesh' ih cherez Berlin?

Davaj ja vse-taki zaedu k tebe v Chemnitz.
Ja ne otnimu u teb'a mnogo vremeni, esli hochesh', prosto zaidu
k tebe nenadolgo.
Ja ponimaju, chto tebe eto nelepym kazhets'a i nenuzhnym,
no mne tak budet gorazdo spokojnee. 

*/


TableFunction::TableFunction(int dim,realtype *arrx,realtype *arry){
 xx=arrx;
 yy=arry;
 n=dim;
 if(n<=0)fatal_error("TableFunction: Invalid array dimension: %d\n",n);
 alloc=0;
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 ftype=XYDATA;
 btype=OUT_ZERO;
}

TableFunction::TableFunction(char *file,int col1=-1,int col2=-1){
 xx=yy=NULL;
 n=ReadTable(file,xx,yy,col1,col2);
 if(n<=0)fatal_error("Can not initialize TableFunction from file %s\n",file);
 if(xx==NULL){
  ftype=YDATA;
  alloc=0x1;
 }
 else{
  ftype=XYDATA;
  alloc=0x3;
 }
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 btype=OUT_ZERO;
}

TableFunction::TableFunction(FILE *f, int col1=-1, int col2=-1){
 xx=yy=NULL;
 n=ReadTableDirect(f,xx,yy,col1,col2,0);
 if(n<=0)fatal_error("Can not initialize TableFunction from file\n");
 if(xx==NULL){
  ftype=YDATA;
  alloc=0x1;
 }
 else{
  ftype=XYDATA;
  alloc=0x3;
 }
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 btype=OUT_LAST;
}


TableFunction& TableFunction::operator<<(const TableFunction &F){
  if((alloc&0x01))delete[] yy;
  if((alloc&0x02))delete[] xx;
  
  memcpy((void*)this,(void*)&F,sizeof(TableFunction));
  
  if(ftype==XYDATA)alloc=0x03;
  else alloc=0x01;
  
  int i;
  if(alloc&0x01){
    yy = new realtype[n];
    if(!yy)serror("TableFunction: MAE.\n");
    for(i=0;i<n;i++)yy[i]=F.yy[i];
  }
  if(alloc&0x02){
     xx = new realtype[n];
     if(!xx)serror("TableFunction: MAE.\n");
     for(i=0;i<n;i++)xx[i]=F.xx[i];
  }
  return *this;
}



void TableFunction::write(char *file, realtype yscale, int (*mapfunc)(realtype *x, realtype *y)){
 int i;
 realtype x=0, y;
 FILE *f=Err_fopen(file,"wt");
 for(i=0;i<n;i++){
   switch(ftype){
   case XYDATA: 
     x=xx[i];
     break;
   case YDATA:
     x=xstart;
     if(n>1)x+=(-xstart+xend)*i/(n-1);
   }
   y=yy[i]*yscale;
   if(mapfunc){
     if(mapfunc(&x,&y)<=0)continue;
   }
   fprintf(f,"" REALFRM " " REALFRM "\n",x,y);
 }
 fclose(f);
}
     

realtype TableFunction::operator()(realtype x){
 
  //printf("" REALFRM " " REALFRM " " REALFRM "\n",x,xend,xstart);
 if(x<xstart){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[0];
  }
 }
 if(x>xend){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[n-1];
  }
 }
 

 realtype d=0;
 int i;

 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;
 if(i>=n-1)return yy[n-1];

 //return yy[0]+x;

 if(ftype==YDATA){
  d=d-(realtype)i;
  d=*(yy+i)+(*(yy+i+1)-*(yy+i))*d;
  return d;
 }
 
 i++;
 while(xx[i]<x)i++;
 while(xx[i-1]>x)i--;

 d=xx[i]-xx[i-1];
 if(d!=0.)d=(x-xx[i-1])/d;
 return yy[i-1]+(yy[i]-yy[i-1])*d;
}

realtype TableFunction::der(realtype x){
 
  //printf("" REALFRM " " REALFRM " " REALFRM "\n",x,xend,xstart);
 if(x<xstart)
   return 0.;
 if(x>xend)
   return 0.;
   
 realtype d=0;
 int i;

 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;
 if(i>=n-1)return 0.;

 //return yy[0]+x;

 if(ftype==YDATA){
   return yy[i+1]-yy[i];
 }
 
 i++;
 while(xx[i]<x)i++;
 while(xx[i-1]>x)i--;
 d=xx[i]-xx[i-1];
 if(d!=0.)
   return (yy[i]-yy[i-1])/d;
 return 0.;
}

TableFunction& TableFunction::operator*=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];

    if(x>=F.xstart && x<= F.xend){
      yy[i]*=F(x);
    }
  }
  return *this;
}


TableFunction& TableFunction::operator/=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];

    if(x>=F.xstart && x<= F.xend){
      realtype del=F(x);
      if(fabs(del)<1e-32){
        if(fabs(yy[i])>1e-32)yy[i]=1e32;
        else yy[i]=0;
      }
      else yy[i]/=del;
    }
  }
  return *this;
}

TableFunction& TableFunction::operator*=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]*=F(x);
  }
  return *this;
}

TableFunction& TableFunction::operator+=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]+=F(x);
  }
  return *this;
}



TableFunction& TableFunction::operator+=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    if(x>=F.xstart && x<= F.xend){
      yy[i]+=F(x);
    }
  }
  return *this;
}



TableFunction& TableFunction::operator+=(realtype c){
 
  int i;
  for(i=0;i<n;i++){
    yy[i]+=c;
  
  }
  return *this;
}

TableFunction& TableFunction::operator-=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]-=F(x);
  }
  return *this;
}



TableFunction& TableFunction::operator-=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    if(x>=F.xstart && x<= F.xend){
      yy[i]-=F(x);
    }
  }
  return *this;
}



TableFunction& TableFunction::operation(TableFunction &F,
					realtype fop(realtype,realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    if(x>=F.xstart && x<= F.xend){
      yy[i]=fop((*this)(x),F(x));
    }
  }
  return *this;
}





TableFunction& TableFunction::operator-=(realtype c){

  int i;
  for(i=0;i<n;i++){
    yy[i]-=c;

  }
  return *this;
}



TableFunction& TableFunction::operator*=(realtype c){

  int i;
  for(i=0;i<n;i++){
    yy[i]*=c;

  }
  return *this;
}


realtype TableFunction::Dx(realtype x){

 if(ftype==YDATA){
   if(n<=1)return 0.;
   else return (xend-xstart)/(n-1);
 }

 if(x<xstart || x>=xend || n<2){
   return 0.;
 }

 realtype d=0;
 int i;
 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;

 i++;
 while(xx[i]<x)i++;
 while(xx[i-1]>x)i--;
 d=xx[i]-xx[i-1];
 if(d!=0.)d=(x-xx[i-1])/d;
 i--;

 realtype dx1, dx2;

 if(i<=0)dx1=xx[1]-xx[0];
 else if(i>=n-1)dx1=xx[n-1]-xx[n-2];
 else dx1=(xx[i+1]-xx[i-1])/2;

 if(i>=n-2)dx2=xx[i+1]-xx[i];
 else dx2=(xx[i+2]-xx[i])/2;

 d=dx1+(dx2-dx1)*d;
 return d;
}

realtype TableFunction::Dy(realtype x){
 if(x<xstart || x>=xend || n<2){
   return 0.;
 }

 realtype d=0;
 int i;

 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;

 if(ftype!=YDATA){
   i++;
   while(xx[i]<x)i++;
   while(xx[i-1]>x)i--;
   d=xx[i]-xx[i-1];
   if(d!=0.)d=(x-xx[i-1])/d;
   i--;
 }
 else{
   d=d-(realtype)i;
 }

 realtype dy1, dy2;

 if(i<=0)dy1=yy[1]-yy[0];
 else if(i>=n-1)dy1=yy[n-1]-yy[n-2];
 else dy1=(yy[i+1]-yy[i-1])/2;

 if(i>=n-2)dy2=yy[i+1]-yy[i];
 else dy2=(yy[i+2]-yy[i])/2;

 d=dy1+(dy2-dy1)*d;
 return d;
}




realtype TableFunction::xscale(realtype x1,realtype x2){
 realtype k/*,x0*/;
 if(x1>x2){k=x1;x1=x2;x2=k;}
 k=(x2-x1)/(xend-xstart);

 int i;
 switch(ftype){
  case XYDATA:
   for(i=0;i<n;i++)xx[i]=x1+k*(xx[i]-xstart);
  case YDATA:
   xstart=x1;
   xend=x2;
 }
 return k;
}

void TableFunction::insert_begin(int steps){
  nsteps=steps;
  stpk=0;
  stpx=0.;
  stpi=0;
  insert_stat=0;
  yold=0.;
}


realtype TableFunction::insert_next(realtype y){
   int dumm;
   realtype d, yr=yold;

   if(stpi==0){ // first step
     yy[stpi++]=y;
     stpx=((realtype)(nsteps-1))*stpi/(n-1); // where the next insert place is
   }
   else{ // scaling
     while((realtype)stpk>=stpx && stpi<n){
       d=DxScale(1.,stpx,dumm);
       yr=yold*(1-d)+y*d;

       yy[stpi]=yr;

       stpi++;
       stpx=((realtype)(nsteps-1))*stpi/(n-1);// where the next insert place is

     }
   }
   if(stpi>=n){
     insert_stat=0;
   }
   else{
     yold=y;
     stpk++; // where the next data comes to
     insert_stat=1;
   }
   return yr;
}


TableFunction * TableFunction::Fcur=NULL;
void SetFunc(TableFunction *f){ TableFunction::Fcur=f;}

realtype TabFunc(realtype x){
 if(TableFunction::Fcur)return (TableFunction::Fcur)->operator()(x);
 else return 0.;
}

realtype TableFunction::integral(){
  if(n<2)return 0.;
  int i;
  realtype sum=0;
  if(ftype==YDATA){
    realtype scale=(xend-xstart)/(n-1);
    for(i=0;i<n-1;i++)sum+=0.5f*(yy[i]+yy[i+1]);
    sum*=scale;
  }
  else{
    for(i=0;i<n-1;i++)sum+=0.5f*(yy[i]+yy[i+1])*(xx[i+1]-xx[i]);
  }
  return sum;
}

TableFunction& TableFunction::integrate(realtype scale){
  if(n<2)return *this;
  int i;
  realtype sum=0, left=0;
  if(ftype==YDATA){
    scale*=(xend-xstart)/(n-1);
    for(i=0;i<n;i++){
      sum+=0.75*yy[i]; // 3/4
      if(i>0)sum+=left/8.;
      else sum+=yy[i]/8.;
      if(i<n-1)sum+=yy[i+1]/8.;
      else sum+=yy[i]/8.;
      left=yy[i];
      yy[i]=sum*scale;
    }
  }
  else{
    double rx=1, lx=1;
    for(i=0;i<n;i++){
      if(i<n-1)rx=0.5*(xx[i+1]-xx[i]);
      if(i==0)lx=rx;
      sum+=0.75*yy[i]*(rx+lx); // 3/4
      if(i>0)sum+=left*lx/8.;
      else sum+=yy[i]*lx/8.;
      if(i<n-1)sum+=yy[i+1]*rx/8.;
      else sum+=yy[i]*rx/8.;
      left=yy[i];
      yy[i]=sum*scale;
      lx=rx;
    }
  }
  return *this;
}

realtype TableFunction::integral_reaches(realtype val){
  int i;
  realtype sum=0;
  if(ftype==YDATA){
    val/=(xend-xstart);
    for(i=0;i<n-1;i++){
      sum+=0.5f*(yy[i]+yy[i+1]);
      if(sum>=val)break;
    }
    if(i>=n-1)return xend;

    realtype sum0=sum-0.5f*(yy[i]+yy[i+1]);
    realtype dy=yy[i+1]-yy[i];

    realtype det=yy[i]*yy[i]+2*dy*(val-sum0);
    if(fabs(dy)<1e-32 || det <0)return xstart+(xend-xstart)*i/(n-1);
    realtype dx=(realtype)(sqrt(det)-yy[i])/dy;
    return xstart+(xend-xstart)*i/(n-1)+dx;
  }
  else{
    for(i=0;i<n-1;i++){
      sum+=0.5f*(yy[i]+yy[i+1])*(xx[i+1]-xx[i]);
      if(sum>=val)break;
    }
    if(i>=n-1)return xend;
    realtype sum0=sum-0.5f*(yy[i]+yy[i+1]);
    realtype dy=yy[i+1]-yy[i];

    realtype det=yy[i]*yy[i]+2*dy*(val-sum0);
    if(fabs(dy)<1e-32 || det <0)return xx[i];
    realtype dx=(realtype)(sqrt(det)-yy[i])/dy;
    return xx[i]+dx;
  }
}



void SetInOrder(unsigned *order,unsigned n, void *array, size_t size){
 unsigned i,k;
 void *tmp=malloc(size);
 if(!tmp)fatal_error("SetInOrder: Memory allocation error.\n");
 char *a1=(char *)array;

 for(i=0;i<n;i++){
  if(order[i]!=i){
   k=order[i];
   while(k<i)k=order[k];
   if(k!=i){
    memcpy(tmp,a1+size*i,size); // swapping
    memcpy(a1+size*i,a1+size*k,size);
    memcpy(a1+size*k,tmp,size);
   }
  }
 }
 free(tmp);
}


realtype *CmpArray;
int CmpAction=1;

int CompIndexed(const void *c1,const void *c2){
 unsigned i1,i2;
 i1=*(int*)c1;
 i2=*(int*)c2;
 if(CmpArray[i1]<CmpArray[i2])return -CmpAction;
 if(CmpArray[i1]>CmpArray[i2])return  CmpAction;
 return 0;
}


void TableFunction::sort(){
 CmpArray=xx;
 unsigned *order=new unsigned[n];
 if(!order)fatal_error("TableFunction:sort: memory allocation error.\n");
 unsigned i;
 for(i=0;i<(unsigned)n;i++)order[i]=i;
 qsort(order,n,sizeof(int),CompIndexed);
 SetInOrder(order,n,xx,sizeof(realtype));
 SetInOrder(order,n,yy,sizeof(realtype));
 for(i=1;i<(unsigned)n;i++){
  if(fabs(xx[i-1]-xx[i])<1e-32){
   memmove(xx+i,xx+i+1,sizeof(int)*(n-(i+1)));
   n--;
  }
 }
 delete [] order;
}


SplineFunction::SplineFunction(int dim,realtype *arr)
 :TableFunction(dim,arr){
 allocspl();
 xx=new realtype[n];
 if(!xx)fatal_error("SplineFunction: Memory allocation error.\n");
 for(int i=0;i<n;i++){
  xx[i]=((realtype)i)/(n-1);
 }
 makespl();
 delete [] xx;
 xx=NULL;
}

SplineFunction::SplineFunction(int dim,realtype *arrx,realtype *arry)
 :TableFunction(dim,arrx,arry){
 allocspl();
 makespl();
}

SplineFunction::SplineFunction(char *file):TableFunction(file){
 allocspl();
 makespl();
}

void SplineFunction::allocspl(){
 aa=new realtype[n-1];
 bb=new realtype[n-1];
 cc=new realtype[n-1];
 if(!aa || !bb || !cc)fatal_error("SplineFunction:allocspl: Memory allocation error.\n");
 alloc|=0x4;
}

// 3-point,initial  conditions: x0*b[0]+x1*c[0]=f[0]
void Passage(int n,realtype *a,realtype *b,realtype *c,realtype *f,realtype *u,realtype *w_arr){
 int i,res=0;
 realtype *l,*k,num;
 if(!w_arr){
  w_arr=new realtype[2*n];
  if(!w_arr)fatal_error("Passage: memory allocation error!\n");
  res=1;
 }
 l=w_arr;
 k=w_arr+n;
 num=1.f/b[0];
 l[0]=-c[0]*num;
 k[0]=f[0]*num;
 for(i=1;i<n;i++){
  num=(b[i]+a[i]*l[i-1]);
  if(fabs(num)<1e-10){
   l[i]=0;
   k[i]=0;
  }
  else{
   num=1/num;
   l[i]=-c[i]*num;
   k[i]=(f[i]-a[i]*k[i-1])*num;
  }
 }
 u[n-1]=k[n-1];
 for(i=n-2;i>=0;i--){
  u[i]=l[i]*u[i+1]+k[i];
 }
 if(res)delete [] w_arr;

}

void SplineFunction::makespl(/*realtype left_der,realtype right_der*/){
 int i;
 realtype d2=0,d3=0,d4=0,dy2=0,dy3=0,dy4=0,k;
 realtype *a_arr[4];

 for(i=0;i<4;i++){
  a_arr[i]= new realtype[n];
  if(!a_arr[i])fatal_error("SplineFunction:makespl: memory allocation error.\n");
 }

 for(i=0;i<n-1;i++){
 /* setting passage coefficients */
  if(i>0){
   d2=d3;
   dy2=dy3;
   d3=d4;
   dy3=dy4;
  }
  if(i<n-2){
   d4=xx[i+2]-xx[i+1];
   dy4=(yy[i+2]-yy[i+1])/d4;
  }

  if(i==0){
   d3=xx[1]-xx[0];
   dy3=(yy[1]-yy[0])/d3;
   a_arr[0][i]=0;
   a_arr[1][i]=2*d3*d4+d3*d3;
   a_arr[2][i]=d4*d4;
   a_arr[3][i]=dy4-2*dy3-dy3*d4/d3;
   continue;
  }

  k=1/(d2+d3);

  if(i==n-2){
   a_arr[0][i]=d3*d2*d2*k;
   a_arr[1][i]=d3*d3*(3-2*k*d3-d2*k);
   a_arr[2][i]=0;
   a_arr[3][i]=(2*d3*(dy2-dy3)-dy2*d3-dy3*d2)*k;
   continue;
  }

  a_arr[0][i]=d2*d2*(d4+d3)*k;
  a_arr[1][i]=3*d3*(d4+d3)-d3*d3*(d4+2*d3)*k-d2*d3*d3*k;
  a_arr[2][i]=d4*d4;
  a_arr[3][i]=((dy2-dy3)*(d4+2*d3)-dy2*d3-dy3*d2)*k+dy4;
 }

 Passage(n-1,a_arr[0],a_arr[1],a_arr[2],a_arr[3],aa,NULL);

 d3=xx[1]-xx[0];
 dy3=(yy[1]-yy[0])/d3;
 bb[0]=dy3/d3 - aa[0]*d3;
 cc[0]=0;

 for(i=1;i<n-1;i++){
  d2=d3;
  dy2=dy3;
  d3=xx[i+1]-xx[i];
  dy3=(yy[i+1]-yy[i])/d3;

  k=1.f/(d2+d3);

  bb[i]=( aa[i-1]*d2*d2-aa[i]*d3*d3+dy3-dy2)*k;
  cc[i]=dy3-bb[i]*d3-aa[i]*d3*d3;
 }

 for(i=0;i<4;i++)delete a_arr[i];
}


realtype SplineFunction::operator()(realtype x){
 if(x<xstart){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[0];
  }
 }
 if(x>xend){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[n-1];
  }
 }
 realtype d;
 int i;
 if(ftype==YDATA){
  d=(x-xstart)/(xend-xstart);
  i=(int)(d*(n-1));
  d=d-(realtype)i/(n-1);
  i++;
 }
 else{
  i=0;
  while(xx[i]<=x)i++;
  if(i>=n)return yy[n-1];
  d=(x-xx[i-1]);
 }
 return ((aa[i-1]*d+bb[i-1])*d+cc[i-1])*d+yy[i-1];
}

realtype SplineFunction::der(realtype x){
 if(x<xstart)
   return 0.;
 if(x>xend) 
   return 0.;
 realtype d;
 int i;
 realtype k=1.;
 if(ftype==YDATA){
   k=(xend-xstart);
  d=(x-xstart)/k;
  i=(int)(d*(n-1));
  d=d-(realtype)i/(n-1);
  i++;
 }
 else{
  i=0;
  while(xx[i]<=x)i++;
  if(i>=n)return yy[n-1];
  d=(x-xx[i-1]);
 }
 return k*((3*aa[i-1]*d+2*bb[i-1])*d+cc[i-1]);
}


realtype SplineFunction::xscale(realtype x1,realtype x2){
 realtype k0=TableFunction::xscale(x1,x2);
 realtype k=1.f/k0;
 int i;
 for(i=0;i<n-1;i++)cc[i]*=k;
 k*=k;
 for(i=0;i<n-1;i++)bb[i]*=k;
 k*=k;
 for(i=0;i<n-1;i++)aa[i]*=k;
 return k0;
}


void SymmMatr::Set(double val){
  int i, n=size*(size+1)/2;
  for(i=0;i<n;i++)arr[i]=val;
}

void SymmMatr::SetDiag(double val){
  unsigned int i;
  for(i=0;i<size;i++)(*this)(i,i)=val;
}

//change_file_strings(char *file,int n, long *pos, char *string[]

# ifdef __BCC
//# pragma argsused
# endif
int no_output(const char *format, ...){
 return 0;
}


int com_div(int a, int b){ // maximal common divider
  a=abs(a);
  b=abs(b);
  int ma=fmax(a,b);
  int mi=fmin(a,b);
  if(mi==0)return 1;
  int rst;
  do{
    rst=ma%mi;
    if(rst==0)return mi;
    rst=ma/mi;
    ma=mi;
    mi=rst;
  }while(mi!=1);

  return mi;
}
 

long com_divl(long a, long b){ // maximal common divider
  a=abs(a);
  b=abs(b);
  long ma=fmax(a,b);
  long mi=fmin(a,b);
  if(mi==0)return 1;
  long rst;
  do{
    rst=ma%mi;
    if(rst==0)return mi;
    rst=ma/mi;
    ma=mi;
    mi=rst;
  }while(mi!=1);

  return mi;
}


long com_divl(int n, long *a){
  int i;
  if(n<=1)return a[0];
  for(i=0;i<n-1;i++){
    a[i]=com_divl(a[i],a[i+1]);
  }
  return com_divl(n-1,a);
}


int cmp_realtypec(const void *c1,const void *c2){
  realtype f1=*(realtype *)c1;
  realtype f2=*(realtype *)c2;
  if(f1<f2)return -1;
  if(f1>f2)return 1;
  return 0;
}



int rand_init=3;

double random1(void){
 static
# ifdef FU
 int
# else
 long
# endif
 i=rand_init;

 i=i*331804469l;
 double R=i*.2328309e-9;
 return R;
}


int get_matching(const char *str,const char *br_pair, const char *skip_cont, const char *quotes, int from_back){
  from_back= from_back ? 1: 0;
  char opening=br_pair[0+from_back], closing=br_pair[1-from_back];
  int ns=0, nq=0;
  if(skip_cont)ns=strlen(skip_cont)/2;
  if(quotes)nq=strlen(quotes);
  int c=1, i=0, s=0, js=-1, q=0, jq=-1;
  int nstr=strlen(str);
  if(from_back)
    str+=nstr-1;
  for(i=0;i<nstr;i++){
    if(!s && *str==opening)c++;
    if(quotes){
      if(!q){ // loking for a match
        for(jq=nq-1;jq>=0;jq--){
          if(*str==quotes[jq]){
            q++;
            break;
          }
        }
      }
      else if(jq>=0){
        if(*str==quotes[jq])q=0;
      }
    }
    if(!q && skip_cont){
      if(!s){ // loking for a match
        for(js=ns-1;js>=0;js--){
          if(*str==skip_cont[2*js+from_back]){
            s++;
            break;
          }
        }
      }
      else if(js>=0){
        if(*str==skip_cont[2*js+from_back])s++;
        if(*str==skip_cont[2*js+1-from_back]){
          if(skip_cont[2*js]==skip_cont[2*js+1])s=0;
          else s--;
        }
      }
    }
    if(!q && !s && *str==closing){
      if(opening==closing)c=0;
      else c--;
    }
    if(c==0)break;
    //i++;
    //str++; 
    if(from_back)str--;
    else str++;
  }
  if(c)return -1;
  return (from_back? nstr-i-1 : i);
}

/*# ifdef UNIX

int random(int num){
 return (int)(num*(0.5+random1()));
}

# endif
*/
# if FU

# ifdef __cplusplus
extern "C" {
#endif

//int vsscanf(const char * str,const char *format,va_list ap){
//  int k=0,i;
//  char *p[15];
//  for(i=0;format[i];i++){
//    if(format[i]=='%'){
//     if(format[i+1]=='%')i++;
//     else k++;
//    }
//  }
//  if(k>13)fatal_error("Ugly vsscanf: can't scan more than 13 fields.\n");
//
//  //va_list ap;
//  //va_start(ap,format);
//  for(i=0;i<k;i++){
//    p[i]=va_arg(ap,char*);
//  }
//
//  switch(k){
//  case 1:
//    return sscanf(str,format,p[0]);
//  case 2:
//    return sscanf(str,format,p[0],p[1]);
//  case 3:
//    return sscanf(str,format,p[0],p[1],p[2]);
//  case 4:
//    return sscanf(str,format,p[0],p[1],p[2],p[3]);
//  case 5:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4]);
//  case 6:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5]);
//  case 7:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
//  case 8:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
//  case 9:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8]);
//  case 10:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9]);
//  case 11:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10]);
//  case 12:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11]);
//  case 13:
//    return sscanf(str,format,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12]);
//  }
//  return 0;
//}

# ifdef __cplusplus
} //extern "C" 
#endif

# endif


# endif









