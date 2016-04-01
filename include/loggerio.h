/*s***************************************************************************
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.4 $
 *   $Date: 2012/06/29 10:50:12 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/loggerio.h,v 1.4 2012/06/29 10:50:12 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/loggerio.h,v $
$Revision: 1.4 $
$Author: valuev $
$Date: 2012/06/29 10:50:12 $
*/
/*s****************************************************************************
 * $Log: loggerio.h,v $
 * Revision 1.4  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.19  2011/09/27 12:36:59  biaks
 * some bags fixed for mingw config
 *
 * Revision 1.18  2011/03/10 02:06:19  lesha
 * em part is moved out from detectors
 *
 * Revision 1.17  2010/02/05 13:56:24  valuev
 * Tuned the Doxyfile for ivutils and added documentation for Plane_3
 *
 * Revision 1.16  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.15  2008/05/06 11:01:03  lesha
 * temporary saving is included
 *
 * Revision 1.14  2008/04/29 01:23:41  lesha
 * Flush is added
 *
 * Revision 1.13  2008/04/23 08:07:16  lesha
 * fflush is commented
 *
 * Revision 1.12  2008/04/08 23:37:39  lesha
 * emFluxSpectra is added
 *
 * Revision 1.11  2008/01/22 10:15:02  lesha
 * unnecessary include files are removed
 *
 * Revision 1.10  2008/01/07 22:29:45  lesha
 * flag is_header is added
 *
 * Revision 1.9  2007/06/22 12:10:28  valuev
 * added file open check
 *
 * Revision 1.8  2007/04/30 19:29:15  lesha
 * _FILE_OFFSET_BITS is commented
 *
 * Revision 1.7  2007/03/25 03:24:56  lesha
 * _FILE_OFFSET_BITS is included
 *
 * Revision 1.6  2007/03/22 14:58:48  lesha
 * bfseek is added
 *
 * Revision 1.5  2007/03/21 01:56:20  lesha
 * long -> long long
 *
 * Revision 1.4  2007/03/20 17:20:38  lesha
 * error with long type in fseek is fixed
 *
 * Revision 1.3  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/
# ifndef LOGGERIO_H
# define LOGGERIO_H

//# ifdef UNIX
//# define _FILE_OFFSET_BITS 64 // for record files > 2 Gb
//# endif

# include <stdio.h>
# include <string>
//#include <sys/stat.h>

using namespace std;

/*e @file loggerio.h
    @brief Classes for storing record/logger sequences in files. 
  
**/

/// class for storing slices of data_t values d:
/// i            <-- sequence length (size_t)
/// i            <-- number of slices (size_t)
/// d d d d d d  <-- slice0
/// d d d d d d  <-- slice1, etc.
/// in a file
/// the slices are of equal length
/// error codes for all functions:
///              -1 out of range
///              -2 invalid file state (not open/ already open)
///              -3 invalid file (can't open)
///              -4 wrong file format
class SeqRecord {
  /// length of unit data in bytes
  size_t data_size;
  /// length of the data sequence
  size_t seqlen;
  /// current file
  FILE *file;
  /// the maximum number left for the last slice
  size_t mnlcs;
  /// the current number of slices
  size_t numslices;
  /// presence of header at the begining of the file
  bool if_header;
  /// file name
  string fname;
public:
  /// constructor
  /// parameter dlen is the length of each slice
  SeqRecord(size_t ds=sizeof(double),size_t dlen=1,bool if_header_=false):
    data_size(ds),seqlen(dlen),file(NULL),if_header(if_header_){}

  /// returns -2 if file is open: can't change existing file
  int SetDataSize(size_t ds) {
    if(file)return -2;
    data_size=ds;
    return 1;
  }

  int SetSeqLen(size_t dlen) {
    if(file)return -2;
    seqlen=dlen;
    return 1;
  }

  /// opens a file
  /// if the mode is "w" then replaces the old file if it exists and starts new record
  /// if the mode is "a" then appends the old file if it exists (after checking format ) or starts new record if not
  /// if the mode is "r" analyses the given file (checks it format) and can perform GetData()
  /// if the mode is "ar" can do both reading and appending
  /// other combinations of r and w are not supported
  /// @return the current slice number, -2 if file is already open

  int OpenRecord(const string &name="record.dat", const char *mode="w", const int shift=0) {
    if(file)return -2; // already open
    file=fopen(name.c_str(),"r");
    bool is_file=(file!=NULL);
    if (is_file)
      fclose(file);
//    struct _stat buf;
//    bool is_file=(!_stat(name.c_str(), &buf));
    char real_mode[10];

    if (*mode=='r') {
      if (is_file) {
        strcpy(real_mode, "rb");
      }
      else {
        return -3;
      }
    }
    else if (*mode=='w') {
      strcpy(real_mode, "w+b");
      is_file=false;
    }
    else if (*mode=='a') {
      if (is_file) {
        strcpy(real_mode, "r+b");
      }
      else {
        strcpy(real_mode, "w+b");
      }
    }
    else
      return -3;
    
    file=fopen(name.c_str(),real_mode);
    if(!file)
      return -3;

    mnlcs=seqlen;

    if (is_file) {
      if (if_header) {
        size_t seqlen1;
        int res=(int)fread(&seqlen1,sizeof(size_t), 1, file);
        res+=(int)fread(&numslices, sizeof(size_t), 1, file);
        if(res!=2 || seqlen1!=seqlen /* fp<2*sizeof(size_t)+sizeof(data_t)*seqlen*numslices*/){
          // wrong format
          CloseRecord();
          return -4; 
        }
      }
      else {
        if(shift){
          numslices=shift;
          long long pos=(long long)(numslices)*(long long)(seqlen)*(long long)(data_size);
//# ifdef UNIX
#ifndef _MSC_VER // значит GCC компилятор
          fseeko64(file, pos, SEEK_SET);
# else
         _fseeki64(file, pos, SEEK_SET);
# endif
        }
        else{
          // checking the size
//# ifdef UNIX
#ifndef _MSC_VER // значит GCC компилятор
          fseeko64(file, 0, SEEK_END);
          long long fsize=ftello64(file);
# else
          _fseeki64(file, 0, SEEK_END);
          long long fsize=_ftelli64(file);
# endif
          numslices=(size_t)(fsize/(data_size*seqlen));
        }
      }
    }
    else {
      numslices=0;
      if (if_header) {
        fwrite(&seqlen,sizeof(size_t),1,file);
        fwrite(&numslices,sizeof(size_t),1,file);
//        fflush(file);
      }
    }
    fname=name;
    return (int)numslices;
  }

  /// appends next arrlen records to a file from rec array
  /// if arrlen exceeds the maximum number left for current slice, appends the correct number of elements only
  /// DATA IS APPENDED BY COPYING FROM MEMORY (using fwrite) <- this behaviour may be changed in the future...
  /// returns the number of elements appended
  int AppendData(const void *rec, size_t arrlen=1) {
    if(!file)return -2;
    size_t count=(mnlcs>arrlen) ? arrlen : mnlcs;
//    long long pos=(if_header ? (long long)(2*sizeof(size_t)) : 0)+((long long)(numslices+1)*(long long)(seqlen)-(long long)(mnlcs))*(long long)(sizeof(data_t));
# ifdef UNIX
//    fseeko(file, pos, SEEK_SET);
# else
//   _fseeki64(file, pos, SEEK_SET);
# endif
    count=fwrite(rec, data_size, count, file);
    mnlcs-=count;
//    fflush(file);
    return (int)count;
//    return fwrite(rec, sizeof(data_t), arrlen, file);
  }

  /// switches to the next slice and fflushes the file
  /// @return the current slice number or -1 in case of error
  int NextSlice() {
    if(!file)return -2;
    if (!mnlcs) {
      numslices++;
      if (if_header) {
        fseek(file, sizeof(size_t), SEEK_SET);
        fwrite(&numslices, sizeof(size_t), 1, file);
//        fflush(file);
      }
      mnlcs=seqlen;
      return (int)numslices;
    }
    else
      return -1;
  }

  int Flush(){
    if(!file)return -2;
    return fflush(file);
  }

  /// closes the file
  /// if the data sequence is not finished, removes the unfinished part from file (OR ?WORSE?: appends default data_t)
  /// @return the current slice number
  int CloseRecord() {
  /// doesn't remove unfinished slice, but it isn't considered next time because number of slices is not incremented
    if(file){
      fclose(file);
      file=NULL;
      return (int)numslices;
    }
    else return -2;
    
  }

  /// gets data from  slice  slicenum starting from beg till end
  /// all counters start from zero
  /// returns the number of elements read
  int GetData(void *rec, size_t slicenum, size_t beg=0, size_t end=0) {
    if(!file)return -2;
    if ((slicenum<numslices)&&(beg<seqlen)) {
      if (end>=seqlen)end=seqlen-1;
      long long pos=(if_header ? (long long)(2*data_size) : 0)+((long long)(slicenum)*(long long)(seqlen)+(long long)(beg))*(long long)(data_size);
//# ifdef UNIX
#ifndef _MSC_VER // значит GCC компилятор
      fseeko64(file, pos, SEEK_SET);
# else
     _fseeki64(file, pos, SEEK_SET);
# endif
//      bfseek(file,pos, SEEK_SET);
      return (int)fread(rec, data_size, end-beg+1, file);
    }
    else
      return -1;
  }

  const char *get_fname(){
    return fname.c_str();
  }
};

# endif
