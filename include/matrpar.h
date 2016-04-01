# ifndef __MATRPAR_H
# define __MATRPAR_H

//e @file matrpar.h
//e @brief procedures to read matrix string-indexed parameters


//e sets delimiter character
char *SetMparDelimiter(char *str);

double MparUnread=1.0003e-32;

char MparDelim[100]="=";

char *SetMparDelimiter(char *str){
  strncpy(MparDelim,str,99);
  return MparDelim;
}

int ReadVectorPar(char *name,int nstr, char *labels, double **values){
  char find[200];
  *values = new double[nstr];
  int i;
  for(i=0;i<nt;i++){
    (*values)[i]=MparUnread;
    sprintf(
  
}



# endif;