/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 1993-2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.9 $
 *   $Date: 2015/02/11 16:04:02 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/src/rparam/rparam.c,v 1.9 2015/02/11 16:04:02 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/src/rparam/rparam.c,v $
$Revision: 1.9 $
$Author: valuev $
$Date: 2015/02/11 16:04:02 $
*/
/*s****************************************************************************
 * $Log: rparam.c,v $
 * Revision 1.9  2015/02/11 16:04:02  valuev
 * restructure, added generalized forces
 *
 * Revision 1.8  2015/01/30 09:08:04  valuev
 * working on sweeps
 *
 * Revision 1.7  2009/04/29 15:25:07  valuev
 * fixed vsscanf return value
 *
 * Revision 1.6  2009/02/10 14:20:37  valuev
 * sync with FDTD project
 *
 * Revision 1.4  2009/01/30 13:54:24  valuev
 * restructured as a library
 *
 * Revision 1.5  2007/07/09 21:27:49  valuev
 * plasma with wave packets
 *
 * Revision 1.3  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "rparam.h"
#include "common.h"

char AcFormat[50];

FILE *cur_file=NULL;
char cur_file_name[250];
char Com_char=';';
char dc[100]=":";
char Scan_end=0;
char Help_char=';';  //Com_char;

char WC_char='*'; // wildcard in format

int stop_on_error=1;

long Get_cur_pos(void){ return ftell(cur_file);}
void Set_cur_pos(long lpos){ fseek(cur_file,lpos,SEEK_SET);}

char Set_delimiter(char c){
 char tmp=dc[0];
 dc[0]=c;
 return tmp;
}

char *Get_delimiter(){
  return dc;
}


FILE *Get_cur_file(void){return cur_file;}

char *Get_cur_filename(void){return cur_file_name;}

void Set_comment_char(char c){Com_char=c;}
void Set_help_char(char c){Help_char=c;}

char  Get_help_char(){ return Help_char;}
char  Get_comment_char(){ return Com_char;}

void Set_wildcard_char(char c){WC_char=c;}
char  Get_wildcard_char(){ return WC_char;}



void Set_stop(int s){ stop_on_error=s;}
int  Get_stop(void){ return stop_on_error;}

void Set_scan_end(char c){ Scan_end=c;}

long get_next_string(char *str, int use_comment){
 long pos;
 int i=0;
 char c;

 if(feof(cur_file))return -1;
 pos=ftell(cur_file);
 do{
  c=(char)fgetc(cur_file);
  if(c=='\n' || c==0x0d || c==0x0a){
# ifndef UNIX   // dos TEXT mode
   if(c==0x0d){
     char c1=(char)fgetc(cur_file);
     if(c1!=0x0a)fputc(c1,cur_file);
   }
# endif
   str[i]=0;
   return pos;
  }
  if(use_comment && c==Com_char)c=0;
  if(Scan_end && c==Scan_end)return -1;
  if(i<1000)str[i++]=c;
 }while(!feof(cur_file));
 str[i-1]=0;
 return pos;
}

# define CASHE_SIZE   1000

struct rparam{
 long pos;
 long helppos;  // additional information begin
 struct rparam *next;
} *Cashe[CASHE_SIZE];

void init_Cashe(void){
 int i;
 for(i=0;i<CASHE_SIZE;i++)Cashe[i]=NULL;
}

int get_code(const char *str){
 int c=0;
 while(*str){
  if(*str==dc[0])return (c&0x00ff)+1;
  if(*str==WC_char)return -((c&0x00ff)+1);
  if(*str!=' ')c+=*str;
  str++;
 }
 return 0; /* no ':' in format string */
}


void free_Cashe(void){
 int i;
 struct rparam *next,*elem;

 for(i=0;i<CASHE_SIZE;i++){
  next=Cashe[i];
  Cashe[i]=NULL;
  while(next!=NULL){
   elem=next->next;
   free(next);
   next=elem;
  }
 }
}

int Open_param_file(const char *file){
 char str[1000];
 long pos, helppos;
 int ind,i/*,j*/;
 struct rparam *elem;

 if(cur_file)fatal_error("Parameter file is already initialized!\n");

 cur_file = fopen(file,"rb");
 if(!cur_file){
  if(stop_on_error==1)fatal_error("Cant open the file '%s'\n", file);
  else if(stop_on_error==-1)msg_error("Cant open the file '%s'\n", file);
  return 0;
 }

 strcpy(cur_file_name,file);
 init_Cashe();

 helppos=-1;
 do{
  i=0;
  pos=get_next_string(str,0);
  if(pos<0)return 1;
  while(str[i]==' ')i++;
  if(helppos<0 && str[i]==Help_char)helppos=pos;
  if(str[i]==Com_char)continue;

  /*
  j=0;
  while(str[i+j]!=':'){
   if(str[i+j]==0)break;
   j++;
  }
  if(!str[i+j])continue;

  ind=str[i]&0x07f;*/
  ind=get_code(str+i);
  //printf("%lx, %s, %d\n",pos,str,ind);
  if(ind==0){
    if(str[i]!=Help_char)helppos=-1;
    continue;
  }
  if(ind<0)ind=-ind;
  ind--;

  /*Err_malloc(elem,sizeof(struct rparam));*/
  elem=(struct rparam *)malloc(sizeof(struct rparam));
  if(!elem)fatal_error("Memory allocation error!\n");
  elem->next=Cashe[ind];
  Cashe[ind]=elem;
  elem->pos=pos;
  elem->helppos=helppos;
  helppos=-1;
 }while(1);
}

void Close_param_file(void){
 free_Cashe();
 if(cur_file)fclose(cur_file);
 cur_file=NULL;
}

long ac_helppos=-1; // position of actual help string in file
long ac_formpos=-1; // position of actual format string in file



char *find_spec(const char *format,char *str,int *form_pos){
 int i,k,ind,l;
 struct rparam *next;
 if(!cur_file)fatal_error("Param File is not initialized.\n");

/* i=0;
 while(format[i]==' ')i++;
 ind=format[i]&0x7f;*/
 ind=get_code(format);
 if(!ind)fatal_error("Incorrect format string %s!\n",format);
 if(ind<0)l=0; //asterik found, complete search
 else l=ind-1;

 ac_helppos=-1;
 ac_formpos=-1;
 do{
  next=Cashe[l];
  while(next!=NULL){
   fseek(cur_file,next->pos,SEEK_SET);
   get_next_string(str,1);
   k=0;
   i=0;
   do{
    while(str[k]==' ')k++;
    while(format[i]==' ')i++;

    if(format[i]==WC_char){
     while(str[k] && str[k]!=dc[0])k++;
     while(format[i] && format[i]!=dc[0])i++;
    }

    if(format[i]!=str[k])break;   /* bad spec */
    if(format[i]==dc[0] || format[i]== 0){
     /* reading */
     if(format[i]==0)fatal_error("Incorrect format specification. There must be %c.\n",dc[0]);
     *form_pos=i+1;
     strncpy(AcFormat,str,k);
     AcFormat[k]='\0';
     ac_helppos=next->helppos;
     ac_formpos=next->pos;
     return str+k+1;
    }

    k++;i++;

   }while(1);
   next=next->next;
  }
  l++;
 }while(ind<0 && l<CASHE_SIZE);
 AcFormat[0]='\0';
 if(stop_on_error>0){
  fatal_error("Can't find specification '%s' in file '%s' \n",format,cur_file_name);
 }
 if(stop_on_error<0){
  msg_error("Can't find specification '%s' in file '%s' \n",format,cur_file_name);
 }
 return NULL;
}

extern int vsscanf(const char *, const char *, va_list);


int rprm(const char *format,va_list ap){
 int i/*,k,ind*/,res;
 //struct rparam *next;
 char arr[1000],*str1,*str;


 if(format[0]=='$'){ /* read in current position */
  do{
   ac_formpos=Get_cur_pos();
   str=arr;
   if(get_next_string(str,1)==-1){
    if(stop_on_error>0)fatal_error("End of file reached in file %s",cur_file_name);
    else if(stop_on_error<0)msg_error("End of file reached in file %s",cur_file_name);
    else return -1;
   }
   while(*str==' ' && *str)str++;
  }while(!*str);

  str1=strstr(format,dc);
  if(str1){
   format=str1+1;
   str1=strstr(str,dc);
   if(str1){
    *str1=0;
    strcpy(AcFormat,str);
    str=str1+1;
   }
   else AcFormat[0]='\0';
  }
  else format=format+1;
  i=0;
  ac_helppos=-1; /* no help assumed for current position */
 }
 else str=find_spec(format,arr,&i);

 if(str){
  if(format[i]=='>'){ /* It is very dangerous to use this option without any
                        args in ... !!! */
   //va_start(ap,format);
   str1=va_arg(ap,char *);
   strcpy(str1,str);
   return 1;
  }
  res=vsscanf(str,format+i,ap);
  if(res<0)
    return 0;
  return res;
 }
 return 0;
}

long Set_position(const char *format){
 char str[1000];
 int i;
 if(!find_spec(format,str,&i))return -1;
 return Get_cur_pos();
}


int Read_param(const char *format,...){
 va_list ap;
 va_start(ap,format);
 return rprm(format,ap);
}

int Read_paramn(int n,const char *format,...){
 va_list ap;
 int st;
 va_start(ap,format);
 st=rprm(format,ap);
 if(st!=n){
  if(stop_on_error>0)fatal_error("Can't read %d parameters '%s' (%d read).\n",n,format,st);
  else if(stop_on_error<0)msg_error("Can't read %d parameters '%s' (%d read).\n",n,format,st);
 }
 return st;
}


int Get_cur_help(char *hbuf,int maxlen){
 long pos;
 char tmp_ch, str[1000];
 int i, strl;

 if(ac_helppos<0){
   *hbuf=0;
   return 0;
 }
 //*hbuf=(char *)malloc(maxlen*sizeof(char));
 //if(!(*hbuf))msg_error("Get_cur_help: MAE\n");
 *hbuf=0;

 //strcpy(*hbuf,"mamochka");
 //return 0;

 pos=Get_cur_pos();
 tmp_ch=Com_char;
 Set_comment_char(-1);
 Set_cur_pos(ac_helppos);

 do{
   get_next_string(str,0);
   i=0;
   while(str[i] && str[i]==' ')i++;
   if(str[i]!=Help_char)break;
   strl=maxlen-strlen(hbuf)-2;
   if(strl<0)break;
   strncat(hbuf,str+i+1,strl);
   strcat(hbuf,"\n");
 }while(1);

 /*
 i=strlen(*hbuf);
 if(i==0){
   free(*hbuf);
   *hbuf=NULL;
 }
 else{
   *hbuf=realloc(*hbuf,(i+1)*sizeof(char));
   if(!(*hbuf))msg_error("Get_cur_help: realloc MAE\n");
 } */
 /* restoring values */
 Set_cur_pos(pos);
 Set_comment_char(tmp_ch);

 return 1;
}



/*
void main(void){
 int g,e,f;
 Open_param_file("tmp.d");
 Read_param("Boba : %d",&g);
 Read_param("Bo*  : %d",&g);
 Read_param(";Boda : %d",&g);
 Read_param("  f  f f  f : %d %d ,%d",&g,&e,&f);
 Close_param_file();
}
*/

