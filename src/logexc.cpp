//# ifdef USE_STDAFX
//# include "stdafx.h"
//# endif

# include "logexc.h"
# include <vector>

using namespace std;

message_logger std_log;

message_logger &message_logger::global(){
  if(!glogp){
    std_log.set_global(true);
  }
  return *glogp;
}

message_logger *message_logger::glogp=NULL;
stdfile_logger default_log("",0,stdout,stderr,vblALLBAD|vblMESS1,vblFATAL,1);

const char *fmt(const char *format,...){
  va_list args;
  va_start(args,format);
  static char buff[1024];
  vsnprintf(buff,1024,format,args);
  va_end(args);
  return buff;
}

std::string fmt_s(const char *format,...)
{
    va_list args;
    va_start(args,format);

    {
        char buff[1024*16];
        if(  (unsigned) vsnprintf(buff,sizeof(buff),format,args)
             <  (unsigned) (sizeof(buff)-1) )
        {
            buff[sizeof(buff)-1]=0;
            va_end(args);
            return buff;
        }
    }

    std::vector<char> buff(1024*32);
    for(; (unsigned) vsnprintf( & buff[0], buff.size(),format,args) >= (unsigned) buff.size()-1 ; )
        buff.resize( buff.size() *2 );

    va_end(args);
    return & buff[0];
}
std::string fmt_va(const char *format, va_list args)
{
    {
        char buff[1024*16];
        if( (unsigned) vsnprintf(buff,sizeof(buff),format,args)
             <  (unsigned) (sizeof(buff)-1) )
        {
            buff[sizeof(buff)-1]=0;
            return buff;
        }
    }

    std::vector<char> buff(1024*32);
    for(; (unsigned) vsnprintf( & buff[0], buff.size(),format,args) >= (unsigned) buff.size()-1 ; )
        buff.resize( buff.size() *2 );

    return & buff[0];
}


//fmt::fmt(const char *format,...){
//  va_list args;
//  va_start(args,format);
//  vsnprintf(buf, sizeof(buf) ,format,args);
//}

