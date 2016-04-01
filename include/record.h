# ifndef RECORD_H
# define RECORD_H


# include <string.h>
# include <stdlib.h>
# include <stdio.h>
# include "common.h"

enum{ FR_TEST, FR_NOWRITE, FR_WRITE, FR_DONE };


class RecFrame{
public:
  long pos;
  char name[50];
  size_t size;
  long fq;
  
  long Nseq;
  long idle;
  long delay;

  int headsize;
  int state[2]; // 0- for writing, 1-- for reading
  long regstp[2];
  long count[2]; // current frame count
  long nstp[2];  //  next valid step for frame

  int spec; // frame specificator

  void sout(char *str){
    sprintf(str,"% 20ld %s %ld %ld %ld %ld %ld &", pos, name,
	    (long)size,fq,Nseq,idle,delay);
  }
  int read(FILE *fp);
};

enum R_STATES{
R_VOID=0,R_INIT=1,R_STEP=2,R_WAITQ=3
};

class Record{
protected:
  long lastst[2];

  long pos0;
  long nsteps;

  char file[250];
  RecFrame *frame;

  long step;
  int state;

  void init_val(){
    nfr=0;
    frame=NULL;
    step=-1;
    nsteps=0;
    state=R_VOID;
    recfp=NULL;
    curpos=0;
    lastst[0]=-1;
    lastst[1]=-1;
  }
  int init_file();
  int open_file();
  int wake_FP();
  void sleep_FP();
  int data_state(long stp,int mode);
  int which_data(int mode); // 0-write, 1-read
  virtual void write_stepmark(long stm);
  long read_stepmark();
  long step_position(long stp);
  int RegisterFrame(RecFrame *frm);

  int nfr;
  FILE *recfp;
  virtual void write_charhead();
  long curpos;

public:
  Record(){
    init_val();
    strcpy(file,"tst.rec");
  }
  Record(char *filename){
    init_val();
    strncpy(file,filename,250);
  }
  int RegisterFrame(char *fn,size_t fsz,long frq,long fdelay=0,
			  long nseq=0,long fidle=-1, int spec_=0);


  virtual ~Record(){
    //if(recfp)fclose(recfp);
    sleep_FP();
    free(frame);
  }


  int SetFileName(char *filename){
    if(state<=R_INIT){
      strncpy(file,filename,250);
      return 1;
    }
    else{
      msg_error("Record: SetFileName(%s): file(%s) already being recorded!\n",
		filename,file);
      return 0;
    }
  }

  virtual int Open(const char *filename);
  virtual int Clear(){
    if(state!=R_VOID){
      sleep_FP();
      free(frame);
      init_val();
      return 1;
    }
    return 0;  
  }
  virtual long Step();
  virtual int Query();
  virtual int Fill(void *framebuff);
  virtual int Get(long stp, int frm, void *buff);
  long getNSteps(){
    return nsteps;
  }
  int getNFrames(){
    return nfr;
  }
  int FrameWithSpec(int spc, int fstart=0);
  int FrameWith(char *str);
  long FrameState(long stp, int frm);
  long FrameFreq(int frm);
  long NRegSteps(long stp, int frm); // tells how many non void steps
  // are left before frame idle period
  // returns 0 if stp is idle, -1 if there is idle time for frame

  // functions to avoid calling time_consuming data_state in Get
  int FixReadStep(long stp);
  int NextReadStep();

};




# endif











