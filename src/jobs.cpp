# include <string.h>
# include "jobs.h"

int JobID::idcount=0;

int JobID::Message(char *msg){
  strncpy(message,msg,249);
  return 0;
}
