# include<string.h>
# include<stdio.h>
# include<stdlib.h> 
# include<time.h> 

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

# ifndef UNIX
# include <io.h>
# endif

# include "record.h"
# include "common.h"



int RecFrame::read(FILE *fp){
  char str[500];
  int i=0;
  int stat=0;
  while(!feof(fp) && i<500){
    str[i]=(char)fgetc(fp);
    if(str[i]=='&'){
      str[i]=0;
      stat=1;
      break;
    }
    i++;
  }
  //printf("rfr (%d): %s\n",i,str);
  if(!stat)return 0;
  headsize=(i+1)*sizeof(char);
  //printf("hs: %d\n",headsize);
  long sz;
  stat=sscanf(str,"%ld %s %ld %ld %ld %ld %ld", &pos, name,
	    &sz,&fq,&Nseq,&idle,&delay);
  if(stat!=7)return 0;
  size=(size_t)sz;
  return 1;
}

int Record::RegisterFrame(char *fn,size_t fsz,long frq,long fdelay,
			  long nseq,long fidle, int spec_){
  if(state>R_INIT){
    // insert here late registration
    msg_error("Record:RegisterFrame: can't register while\n"
	      "recording file (%s)!\n",file);
    return -1;
  }

  if(nfr==0){
    frame=(RecFrame *)malloc(sizeof(RecFrame));
    state=R_INIT;
  }
  else{
    frame=(RecFrame *)realloc(frame,(nfr+1)*sizeof(RecFrame));
  }
  if(!frame)fatal_error("Record: RegisterFrame: MAE.\n");

  strcpy(frame[nfr].name,fn);
  frame[nfr].size=fsz;
  frame[nfr].fq=frq;
  frame[nfr].Nseq=nseq;
  frame[nfr].idle=fidle;
  frame[nfr].delay=fdelay;

  frame[nfr].spec=spec_;
  nfr++;
  return nfr;
}

int Record::RegisterFrame(RecFrame *frm){
  
  int res=RegisterFrame(frm->name,frm->size,frm->fq,frm->delay,frm->Nseq,frm->idle);
  if(res>0){
    frame[nfr-1].headsize=frm->headsize;
  }
  return res;
}


int Record::open_file(){
  if(recfp!=NULL){
    printf("Record: open: Warning: file was not closed!\n");
    fclose(recfp);
    recfp=NULL;
  }

# ifndef UNIX
  int fhandle=_open(file,O_BINARY|O_RDWR);

  //recfp=fopen(file,"rb");
# else
  int fhandle=open(file,O_RDWR);
# endif
  if(fhandle==-1)recfp=NULL;
# ifndef UNIX  
  else recfp=_fdopen(fhandle,"w+");
# else 
  else recfp=fdopen(fhandle,"w+");
# endif

  if(!recfp){
    perror("Record");
    msg_error("Record: Can't (re)open file '%s'.\n",file);
    return 0;
  }

  if(fseek(recfp,curpos,SEEK_SET)){
    msg_error("Record: fseek in file '%s' failed!\n",file);
    return 0;
  }
  return 1;
}


int Record::init_file(){
  //printf("now\n");


# ifndef UNIX
  /*int handle=open(file,O_BINARY|O_RDWR|O_CREAT);
  if(handle<0){
    perror("Record");
  }
  recfp=fdopen(handle,"w+");*/

  // here existance check !!!
  /*if(!access(file,0)){ // file exists
    recfp=fopen(file,"r+b");
  }
  else recfp=fopen(file,"w+b");*/
  recfp=fopen(file,"w+b");

# else
  recfp=fdopen(open(file,O_RDWR|O_CREAT,S_IRWXU),"w+");
# endif



  if(!recfp){
    perror("Record");
    msg_error("Record: Can't open file '%s' for writing.\n",file);
    return 0;
  }

  write_charhead(); // writes the header
  fprintf(recfp,"@Record_begin:");
  
  pos0=ftell(recfp);
  

  long pos;
  int i;

  // writing all existent frame headers
  // one after one
  // pos points to the next frame header (in file)
  // and to the actual file pos (in memory) !
  // -1 ends the sequence (in file)

  char str[500];
  unsigned int len;

  for(i=0;i<nfr;i++){
    // printf("now\n");
    pos=ftell(recfp);
    frame[i].sout(str);
    len=frame[i].headsize=strlen(str);
    if(i!=nfr-1)frame[i].pos=pos+len;
    else frame[i].pos=-1;
    frame[i].sout(str);
    if(len!=strlen(str)){
      fatal_error("Record: field length error (check sizeof long)!\n");
    }
    fprintf(recfp,"%s",str);
    frame[i].pos=pos;
  }

  //sleep_FP();
  //open_file(); // opening in read-write mode
  curpos=ftell(recfp);
  write_stepmark(0);
  //printf("ok\n");
  return 1;
}


void Record::write_stepmark(long stm){
  nsteps=stm;
  //printf("seek mark\n");
  fseek(recfp,pos0,SEEK_SET);
  //printf("OK");
  fprintf(recfp,"%10ld",stm);
  fseek(recfp,curpos,SEEK_SET);
}

long Record::read_stepmark(){
  long stm;
  fseek(recfp,pos0,SEEK_SET);
  if(fscanf(recfp,"%ld",&stm)!=1){
    msg_error("Record: read_step: cannot read number of steps in record!\n");
    return 0;
  }
  fseek(recfp,pos0+10*sizeof(char),SEEK_SET);
  nsteps=stm;
  return stm;
}



long Record::step_position(long stp){
  if(stp<0)return -1;
  int i;
  long N;
  long pos=pos0+10*sizeof(char);
  stp--;
  //printf("stpos: %ld\n",pos);

  if(stp>=0){
    for(i=0;i<nfr;i++){   
      long t=stp-frame[i].delay;
      if(t<0)continue; // frame not entered yet
      
      long eff_step;
      if(frame[i].Nseq>0){
	eff_step=(frame[i].Nseq-1)*frame[i].fq;
	if(frame[i].idle>=0)eff_step+=frame[i].idle;
	else{
	  if(t>eff_step){
	    pos+=frame[i].Nseq*frame[i].size; // frame entered Nseq times and exited
	    continue;
	  }
	}
	N=t/(eff_step+1);
	pos+=N*frame[i].Nseq*frame[i].size; // Nseq repeated N times already
	eff_step=t%(eff_step+1);
      }
      else eff_step=t;
      
      N=eff_step/frame[i].fq+1;
      //if(eff_step!=0)N++;
      
      if(frame[i].Nseq>0){
	if(N>frame[i].Nseq)N=frame[i].Nseq;
      }
      pos+=N*frame[i].size; // frame repeated N times
      //printf("stpos: %ld\n",pos);
    }
  }

  for(i=0;i<nfr;i++){
    if(frame[i].pos<=pos)pos+=frame[i].headsize;
  }
  //printf("stpos: %ld\n",pos);
  return pos;
}

int Record::data_state(long stp, int mode=0){
  if(stp<0)return -1;
  if(lastst[mode]==stp)return 0;
  //if(stp-lastst[mode]==1){
    //printf("tr ");
  //NextReadStep();
  //}
  // printf("d ");

  int i;
  for(i=0;i<nfr;i++){
    frame[i].state[mode]=FR_NOWRITE;
    frame[i].regstp[mode]=0;
    long t=stp-frame[i].delay;

    frame[i].nstp[mode]=frame[i].delay;
    if(t<0)continue;
    
    long eff_step;
    if(frame[i].Nseq>0){
      eff_step=(frame[i].Nseq-1)*frame[i].fq;
      if(frame[i].idle>=0)eff_step+=frame[i].idle;
      else{
	if(t>eff_step){
	  frame[i].nstp[mode]=-1;
	  continue;
	}
      }
      eff_step=t%(eff_step+1);
      frame[i].nstp[mode]+=t-eff_step;
      frame[i].count[mode]=t/(eff_step+1);
      frame[i].count[mode]*=frame[i].Nseq;
    }
    else{
      eff_step=t;
      frame[i].regstp[mode]=-1;
      frame[i].count[mode]=0;
    }
    
    long rest=eff_step%frame[i].fq;
    long N=eff_step/frame[i].fq+1;
    frame[i].count[mode]+=N;
    
    if(frame[i].Nseq>0){
      long Nrest=frame[i].Nseq-N+1;
     

      if(Nrest<=1){
	frame[i].nstp[mode]+=frame[i].Nseq*frame[i].fq+frame[i].idle;
      }
      else frame[i].nstp[mode]+=eff_step-rest+frame[i].fq;

      if(Nrest<=0){	
	frame[i].regstp[mode]=0;
	continue;
      }
      frame[i].regstp[mode]=Nrest;
      //if(N>frame[i].Nseq)continue;
    }
    else frame[i].nstp[mode]+=eff_step-rest+frame[i].fq;

    if(rest!=0)continue;
    frame[i].state[mode]=FR_WRITE;
    //printf("set wr fr. %d / %d\n",i,nfr);
  }
  lastst[mode]=stp;
  return 0;
}



long Record::NRegSteps(long stp, int frm){
  if(stp>=nsteps || stp<0)return 0; // invalid step
  if(frm>=nfr || frm<0)return 0; // invalid frame

  data_state(stp,1);
  return frame[frm].regstp[1];
}


long Record::FrameFreq(int frm){
  if(frm>=nfr || frm<0)return 0; // invalid frame
  return frame[frm].fq;
}





int Record::FrameWithSpec(int spc, int fstart){
  if(fstart<0)return -1; // invalid frame
  int i;
  for(i=fstart;i<nfr;i++){
    if(frame[i].spec==spc)return i;
  }
  return -1;
}
 


int Record::which_data(int mode=0){
  if(step<0)return -1;
  int i;
  for(i=0;i<nfr;i++){
    if(frame[i].state[mode]==FR_DONE||frame[i].state[mode]==FR_NOWRITE)
      continue;
    if(frame[i].state[mode]==FR_WRITE)return i;
  }
  return -1;
}


long Record::Step(){
  //printf("stp %ld\n",step);
  if(state!=R_INIT && state!=R_STEP){
    msg_error("Record:Step: %ld cannot step:"
	      "no frames initialized or not all frames are filled yet!\n",step);
    return -1;
  }
  int res;
  if(step==-1){
    res=init_file();
    /*
    int i;
    char str[200];
    for(i=0;i<nfr;i++){
      frame[i].sout(str);
      printf("%s\n",str);
    }*/
  }
  else {
    res=wake_FP();
    //res=open_file();
  }
  // printf("init\n");

  if(!res)return -1;

  step++;
  write_stepmark(step);
  //nsteps=step;

  
  //printf("Step pos: %ld  %ld\n",step_position(step),ftell(recfp));

  data_state(step);
  if(which_data()<0){
    state=R_STEP;
    sleep_FP();
    //write_stepmark(step);
  }
  else state=R_WAITQ;
  
  return step;
}

int Record::Query(){
  if(state!=R_WAITQ)return -1; // no frame awaited
  int res=which_data();
  if(res==-1){
    state=R_STEP;
    //write_stepmark(step);
  } 
  return res;
}


int Record::wake_FP(){
  if(recfp==NULL){
    return open_file();
  }
  else{
    if(fseek(recfp,curpos,SEEK_SET)){
      msg_error("Record: fseek in file '%s' failed!\n",file);
      return 0;
    }
    return 1;
  }
}

void Record::sleep_FP(){
  if(recfp)fclose(recfp);
  recfp=NULL;
}


int Record::Fill(void *framebuf){
  //printf("f1\n");
  int fr=which_data();
  if(fr<0){
    msg_error("Record:Fill: %ld - no data awaited now!\n",step);
    return -1;
  }

  //printf("f2\n");
  state=R_WAITQ;
  if(!wake_FP())return -1;
  
  //printf("f3\n");
  //printf("wr: %ld\n",ftell(recfp));
  if(fwrite(framebuf,1,frame[fr].size,recfp)!=frame[fr].size){
    msg_error("Record: write error!\n");
    return -1;
  }
  frame[fr].state[0]=FR_DONE;
  curpos=ftell(recfp);
  
  
  fr=which_data();
  if(fr<0){
    write_stepmark(step+1);
    sleep_FP();
    state=R_STEP;
  }
  return fr;
}


void Record::write_charhead(){
  fprintf(recfp,"//consecutive periodic record\n");
  time_t t;
  time(&t);
  fprintf(recfp,"Creation date: %s",ctime(&t));
}
  


int Record::Open(const char *filename=NULL){
  if(!filename){
    if(state!=R_VOID){
      eprintf("Record:Open: file '%s' already registered"
	      " and being recorded!\n",file);
      return 0;
    }
  }
  else{
    if(state==R_WAITQ){
      eprintf("Record:Open: warning: closing record '%s'"
	      " with unfinished step!\n",file);
    }
    sleep_FP();
    if(frame)free(frame);
    init_val();
    strcpy(file,filename);
  }

  if(!open_file())return -1;  
  //printf("oo\n");
  // initializing record
  fseek(recfp,0,SEEK_SET);

  char str[500]="Record_begin:";
  int stat=0,i;
  while(!feof(recfp)){
    char ch;
    ch=(char)fgetc(recfp);
    //printf("%c",ch);
    if(ch=='@'){
      stat=1;
      for(i=0;str[i];i++){
	if(str[i]!=fgetc(recfp)){
	  stat=0;
	  break;
	}
      }
      break;
    }
  }
  if(!stat){
    msg_error("Record:Open: file '%s' has invalid format!\n",file);
    return -1;
  }
  //printf("oo1\n");
  pos0=ftell(recfp);
  nsteps=read_stepmark();
  step=nsteps-1;

  // reading frames
  long pos;
  RecFrame tmpfr;
  do{
    pos=ftell(recfp);
    if(!tmpfr.read(recfp)){
      msg_error("Record:Open: Can't read frame #%d from file '%s'\n",
		nfr+1,file);
      return -1;
    }
   
    RegisterFrame(&tmpfr);
    frame[nfr-1].pos=pos;

  }while(tmpfr.pos!=-1);

  state=R_STEP;
  //printf("oo2\n");

  // analyzing data size
  pos=step_position(step+1);
  curpos=pos;
  fseek(recfp,0,SEEK_END);
  //printf("sizes: %ld %ld, pos0: %ld\n",ftell(recfp),pos,pos0);
  
  /*
  printf("open check:\n");
  for(i=0;i<nfr;i++){
    frame[i].sout(str);
    printf("%s hs: %d\n",str,frame[i].headsize);
  }*/

  

  if(ftell(recfp)!=pos){
    fseek(recfp,curpos,SEEK_SET);
    return 1; // unfinished steps at the end
  }
  
  return 0; // everything is OK

}


int Record::FixReadStep(long stp){
  return data_state(stp,1);
}


int Record::NextReadStep(){
  long stp=++(lastst[1]);
  int i;
  for(i=0;i<nfr;i++){
    if(stp<frame[i].nstp[1])frame[i].state[1]=FR_NOWRITE;
    else if(frame[i].nstp[1]==stp){
     
      frame[i].state[1]=FR_WRITE;

      // calculating new nstp     
      
      long tnstp=stp+frame[i].fq;

      if(frame[i].Nseq>0){
	long N=frame[i].count[1]%frame[i].Nseq;
	if(N==0){ // nseq ended
	  if(frame[i].idle>=0){
	    tnstp+=frame[i].idle;
	  }
	  else tnstp=-1;
	  frame[i].regstp[1]=-1;
	}
	else frame[i].regstp[1]=frame[i].Nseq-N;
      }
      else frame[i].regstp[1]=0;
      (frame[i].count[1])++;
      frame[i].nstp[1]=tnstp;
    }
  }
  return stp;
}



long Record::FrameState(long stp, int frm){
  if(stp>=nsteps || stp<0)return 0; // invalid step
  if(frm>=nfr || frm<0)return 0; // invalid frame

  data_state(stp,1);
  if(frame[frm].state[1]!=FR_WRITE)return 0; // no such frame on this step
  else return frame[frm].size;
}
  

int Record::Get(long stp, int frm, void *buff){
  //printf("nstp: %d\n",nsteps);
  if(stp>=nsteps || stp<0)return 0; // invalid step
  if(frm>=nfr || frm<0)return 0; // invalid frame
  
  data_state(stp,1);
  if(frame[frm].state[1]!=FR_WRITE){
    //printf("no frame %d\n",frm);
    return 0; // no such frame on this step
  }

  long pos=step_position(stp);
  if(pos<0)return 0; // invalid step
  //printf("Get pos: %ld\n",pos);


  int i;
  for(i=0;i<frm;i++){
    if(frame[i].state[1]==FR_WRITE)pos+=frame[i].size;
  }

  wake_FP();
  if(fseek(recfp,pos,SEEK_SET)){
    fatal_error("Fseek failed!\n");
  }


  if(fread(buff,1,frame[frm].size,recfp)!=frame[frm].size)
    return 0; // file read failure*
  /*fseek(recfp,0x2a3,SEEK_SET);
  char zzz[100];
  fread(zzz,1,100,recfp);*/


  return 1;
}





int Record::FrameWith(char *str){
  int i;
  for(i=0;i<nfr;i++){
    if(strstr(frame[i].name,str))return i;
  }
  return -1;
}



# if 0

int main(){
  Record Rec;

  Rec.SetFileName("test.rec");
  Rec.RegisterFrame("energy",sizeof(float),3,5,10,-1);
  Rec.RegisterFrame("freq",sizeof(float),2);
  //Rec.RegisterFrame("array",10*sizeof(int),3);

  int n=50;
  int i;
  int array[10];
  float E,E1=-10;

  printf("record:\n");
  for(i=0;i<n;i++){
    //printf("step %d\n",i+1);
    array[i%10]=i;
    E=20.*i;

    printf("step %d\n",i);
    Rec.Step();
    
    int fr=Rec.Query();
    //printf("s\n");
    while(fr>=0){
      if(fr==0){
	fr=Rec.Fill(&E1);
	printf("rec E1: %f\n",E1);
	//printf("q1\n");
	continue;
      }
      else{
	fr=Rec.Fill(&E);
	printf("rec E: %f\n",E);
	//printf("q2\n");
	//fr=Rec.Fill(array); 
	continue;
      }
     
    }
  }

  int j;
  printf("now open:\n");
  j=Rec.Open("test.rec");
  printf("done %d\n",j);
  for(i=0;i<n;i++){
    printf("step #%d:\n",i);
    if(Rec.Get(i,1,&E)){
      printf("out E: %f\n",E);
    }
    if(Rec.Get(i,0,&E1)){
      printf("out E1: %f\n",E1);
    }
    /*if(Rec.Get(i,1,array)){
      printf("array: ");
      for(j=0;j<10;j++){
	printf("%d ",array[j]);
      }
      printf("\n");
    }*/
  }
}
# endif
