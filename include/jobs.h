# ifndef __JOBS_H
# define __JOBS_H

class JobID;

typedef int (*pExeFunc)(void *obj, float &progress);
typedef int (*pNotifyFunc)(JobID *job);


class JobID{

public:
  /*int submit(void *object, pJobFunc func);
  int execute(){
    if(!(status&JOB_SUBMITTED))return 0;
    res=job_func(obj,
  }*/

  JobID(){
    id=idcount++;
    source=NULL;
    target=NULL;
    status=0;
    progress=0;
    frac=1;
  }

  void SetNotificator(pNotifyFunc func){
    notf=func;
  }
  void SetTarget(void *ptr){
    target=ptr;
  }
  void *GetTarget(){
    return target;
  }
  float GetProgress(){
    return progress;
  }
  void SetProgress(float p){
    progress=p;
  }
  int CheckTarget(){
    if(!notf)return -1;
    return (*notf)(this);
  }
  float Fraction(float f){
    frac=f;
    if(frac<=0 || frac>1)frac=1;
    return frac;
  }
  float AddToProgress(float dp){
    progress+=dp*frac;
    return progress;
  }
  char *GetMessage(){
    return message;
  }
  int Message(char *msg);

protected:
  pNotifyFunc notf;
  static int idcount;
  int id;
  int status;
  void *source;
  void *target;
  float progress;
  char message[250];
  float frac;
};


/*
class JobHandler{

public:
  int notify(JobID *job){



} */



# endif
