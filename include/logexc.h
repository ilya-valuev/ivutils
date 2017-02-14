# ifndef LOGEXC_H
# define LOGEXC_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdexcept>
#include <string>

# ifndef _WIN32
# include <typeinfo>
# endif


/// this specifies whether to put file/line info in error messages
# ifndef LINFO
# define LINFO 1
# endif

/// verbosity levels
enum vbLEVELS{
  vblNONE   = 0,  ///< completely silent
  vblFATAL  = 0x1,  ///< report fatal errors
  vblOERR   = 0x2, ///< report other errors
  vblERR    = vblFATAL|vblOERR, ///< report all errors
  vblWARN   = 0x4, ///< report warnings
  vblALLBAD = vblWARN|vblERR, ///< report all errors and warnings
  vblMESS1  = 0x8, ///< report messages of level 1 (important)
  vblMESS2  = 0x10, ///< report messages of level 2 (less important)
  vblMESS3  = 0x20, ///< report messages of level 3 (not important)
  vblMESS4  = 0x40, ///< report messages of level 4 (anything possible)
  vblALLMESS = vblMESS1|vblMESS2|vblMESS3|vblMESS4, ///< report all messages
  vblPROGR  = 0x80,           ///< report progress
  vblALL    = 0xffff
};



/// traits structure to deduce exception level and name from exception class
/// by default all exceptions have vblFATAL level
template<class exc_t>
struct log_exception_traits{
  /// exeption level according to the vbLEVELS
  static int level(const exc_t &signal){ return vblFATAL; } 
  /// the string name of exception category
  static std::string name(const exc_t &signal){ return typeid(exc_t).name();}
  /// adds some more explanations to the description
  /// default behaviour: nothing done
  static exc_t add_words(const exc_t &orig, const char *words){
    return orig;
  }
};




/// integer exceptions have the level equal to their value
template<>
struct log_exception_traits<int>{
  /// exeption level according to the vbLEVELS
  static int level(const int &signal){ return signal; } 
  /// the string name of exception category
  static std::string name(const int &signal){ 
    if(signal&vblFATAL)
      return "fatal error";
    else if(signal&vblOERR)
      return "error";
    else if(signal&vblWARN)
      return "warning";
    else
      return "";
    /*
    else if(signal&vblALLMESS)
      return "message";
    else if(signal&vblPROGR)
      return "progress report";
    else
      return "integer exception";*/
  }
  /// default behaviour: nothing done
  static int add_words(const int &orig, const char *words){
    (void)words;
    return orig;
  }
};

/// vbLEVELS exceptions act as integers
template<>
struct log_exception_traits<enum vbLEVELS>{
  static int level(const enum vbLEVELS &signal){ return log_exception_traits<int>::level(signal); } 
  static std::string name(const enum vbLEVELS &signal){ return log_exception_traits<int>::name(signal); } 
  static enum vbLEVELS add_words(const enum vbLEVELS &orig, const char *words){
    (void)words;
    return orig;
  }
};


// format a string
const char *fmt(const char *format,...);
std::string fmt_s(const char *format,...);
std::string fmt_va(const char *format, va_list);

# if USE_WX
///\en Here arguments are wxStrings (not working!)
const char *wxfmt(const char *format,...);

# endif

/// Logger class to control (computational) function behaviour when something requiring user attention has happened.
/// message(signal,errcode, text) is used to either throw an exception or return errorcode
/// At first, the the level of error is determined via log_exception_traits<>::level(signal)
/// For integer (enum) signals the level is the signal itself.
/// Then text is printed, if signal level is listed in output levels or (or in extra outlevels, if they are set) 
/// via log_text() function.
/// If level has vblERR bit, the behaviour is controlled by the flag specified in set_throw(flag):
/// flag=0:   nothing done;
/// flag=1:   calls add_words() for signal and throws signal;
/// flag=2:   throws pair<>(errcode, text);
/// flag=3:   throws errcode.
/// Then, if the level is listed in stop_levels (or in extra stop levels, if they are set), the program is aborted,
/// otherwise errcode is returned;
/// The function set_levels(out_levels,stop_levels) is used to specify bitflags for the levels which
/// require message output or/and program termination. Stop level has effect only when exceptions are not thrown.
/// The function extra_levels(eout_levels,estop_levels) is used to temporarily set the corresponding levels,
/// they are unset (the original levels are restored) by calling extra_levels(0,0).
class message_logger {
  // global message is a friend
 // template<class exc_t>
  //friend int message(const exc_t &signal, int errcode, const char *what, ...);
#ifdef __MSVC__
  template<class exc_t, typename arg1_type>
  friend int message(const exc_t &signal, int errcode, const char *what, arg1_type, ...);
  template<class exc_t>
  friend int message(const exc_t &signal, int errcode, const char *what);
#else
    template<class exc_t>
    friend int message(const exc_t &signal, int errcode, const char *what, ...);
    template<class exc_t>
    friend int message_str(const exc_t &signal, int errcode, const char *what);
#endif
protected:
  std::string descriptor;
  int throw_ex;
  int outlevel, eoutlevel;
  int stoplevel, estoplevel;
  
  static message_logger *glogp;
  /// used to restore the previous global logger
  message_logger *prev, *next;
public:
  
  message_logger(const std::string &descriptor_="", int out_level=vblALLBAD|vblMESS1, 
                 int stop_level=vblFATAL, int throw_exceptions=0, int use_globally=0)
    :descriptor(descriptor_),prev(NULL),next(NULL){
    set_throw(throw_exceptions);
    set_levels(out_level,stop_level);
    extra_levels(0,0);
    if(use_globally){
      set_global(true);
    }
  }
  
  /// returns a reference to global logger
  /// if not set, links with default message_logger
  static message_logger &global();
  
  /// sets/unsets this logger as the global logger
  int set_global(bool set){
    if(set){
      if(prev) // already set
        return -1;
      if(glogp)
        glogp->next=this;
      prev=glogp;
      glogp=this;
    }
    else{
      if(glogp!=this) // was not set as the global
        return -1; 
      glogp=prev;
      if(glogp)
        glogp->next=NULL;
      prev=NULL;
    }
    return 1;
  }
  
  virtual void set_throw(int throw_exceptions){
    throw_ex=throw_exceptions;
  }

  virtual void set_levels(int out_level=vblALLBAD|vblMESS1, int stop_level=vblFATAL){
    outlevel=out_level;
    stoplevel=stop_level;
  }

  /// nonzero extra levels are applied instead of set ones
  virtual void extra_levels(int out_level=vblALLBAD|vblMESS1, int stop_level=vblFATAL){
    eoutlevel=out_level;
    estoplevel=stop_level;
  }
  

  template<class exc_t>
  int message(const exc_t &signal, int errcode, const char *what){
      return str_message(signal, errcode, what);
  }

  template<class exc_t, typename arg1_type>
  int message(const exc_t &signal, int errcode, const char *what, arg1_type arg1, ...){
    (void) arg1;
    va_list args;
    va_start(args,what);
    int ret= va_message(&signal, errcode, what, args);
    va_end(args);
    return ret;
  }

  template<class exc_t>
  int str_message(const exc_t &signal, int errcode, const char *what){
      int level=log_exception_traits<exc_t>::level(signal);
      if(level&(eoutlevel ? eoutlevel : outlevel)){ //needs to print a message
        log_text(level,(const char *)log_exception_traits<exc_t>::name(signal).c_str(), what);
      }
      if(level&vblERR){
        if(throw_ex==1) // throws exc_t exception
          throw log_exception_traits<exc_t>::add_words(signal,what);
        else if(throw_ex==2) // throws pair<>(int,const char*) exception
          throw std::make_pair(errcode,what);
        else if(throw_ex==3) // throws int exception
          throw errcode;
      }
      if(level&(estoplevel ? estoplevel: stoplevel) ){ // needs to stop
        exit(errcode);
      }
      return errcode;
  }

  template<class exc_t>
  int va_message(const exc_t &signal, int errcode, const char *what, va_list args){
    int level=log_exception_traits<exc_t>::level(signal);

    if(level&(eoutlevel ? eoutlevel : outlevel)){ //needs to print a message
      return message(signal, errcode,  fmt_va(what, args).c_str() );
    }
    else {
      return message(signal, errcode, what);
    }
  }

  virtual void log_text(int level, const char *messtype, const char *messtext){
    (void)level;
    if(descriptor!="") // descriptor is used as header
      printf("%s:\n",descriptor.c_str());
    //if(messtype!="" )
    if(messtype!=0 && messtype[0]!=0 )
      printf("%s: ",messtype);
    printf("%s\n",messtext);
  }

  /// checks that the deleted one is not in global logger chain 
  ~message_logger(){
    if(prev){
      prev->next=next;
      if(next)
        next->prev=prev;
    }
    set_global(false);
  }
};

#ifdef __MSVC__
/// global message function
template<class exc_t, typename arg1_type>
int message(const exc_t &signal, int errcode, const char *what, arg1_type arg1, ...){
  if(message_logger::glogp){
    (void) arg1;
    va_list args;
    va_start(args,what);
    return message_logger::glogp->va_message(signal,errcode,what,args);
  }
  else 
    return -1;
}

/// global message function
template<class exc_t>
int message(const exc_t &signal, int errcode, const char *what){
  if(message_logger::glogp){
    return message_logger::glogp->message(signal,errcode,what);
  }
  else
    return -1;
}
#else
template<class exc_t>
int message(const exc_t &signal, int errcode, const char *what, ...){
  if(message_logger::glogp){
    va_list args;
    va_start(args,what);
    int ret=message_logger::glogp->va_message(signal,errcode,what,args);
    va_end(args);
    return ret;
  }
  else
    return -1;
}
#endif

template<class exc_t>
int message_str(const exc_t &signal, int errcode, const char *what){
  if(message_logger::glogp){
    return message_logger::glogp->str_message(signal,errcode,what);
  }
  else
    return -1;
}

/// message logger for which std and error streams may be specified
class stdfile_logger: public message_logger {
protected:
  FILE *fout, *ferr;
public:
  stdfile_logger(const std::string &descriptor_="", int throw_exceptions=0, 
                  FILE *out=stdout, FILE *err=stderr, 
                  int out_level=vblALLBAD|vblMESS1,int stop_level=vblFATAL,
                  int use_globally=0)
                  : message_logger(descriptor_,out_level,stop_level,throw_exceptions,use_globally),fout(NULL), ferr(NULL){
    set_out(out);
    set_err(err);
  }
  virtual void set_out(FILE *out, int close_prev=0){
    if(close_prev && fout && fout!=stderr && fout !=stdout)
      fclose(fout);
    fout=out;
  }
  
  virtual void set_err(FILE *err, int close_prev=0){
    if(close_prev && ferr && ferr!=stderr && ferr !=stdout)
      fclose(ferr);
    ferr=err;
  }

  virtual void log_text(int level, const char *messtype, const char *messtext){
    FILE *f= (level&vblALLBAD ? ferr : fout);
    if(!f)
      return;
    if(descriptor!="") // descriptor is used as header
      fprintf(f,"%s:\n",descriptor.c_str());
    if(std::string(messtype)!="")
      fprintf(f,"%s: ",messtype);
    fprintf(f,"%s\n",messtext);
  }
};


/// macros with common usage
#define LOGFATAL(code,text,lineinfo) ((lineinfo) ? ::message(vblFATAL,(code)," %s at %s:%d",(text),__FILE__,__LINE__) : \
                                   (::message(vblFATAL,(code), (text))) )


#define LOGERR(code,text,lineinfo) ((lineinfo) ? ::message(vblOERR,(code)," %s at %s:%d",(text),__FILE__,__LINE__) : \
                                   (::message(vblOERR,(code), (text))) )


#define LOGERR0(text) LOGERR(0, (text), 1)


#define LOGMSG(cat,text,lineinfo) ((lineinfo) ? ::message((cat),0," %s at %s:%d",(text),__FILE__,__LINE__) : \
                                  (::message((cat),0, (text))) )



# if 0 /// did not work

/// this may be used to inherit exceptions
/// where level and name are defined whithin a class
struct log_exception {
   /// exeption level according to the vbLEVELS
  static int level(const log_exception &signal){ return vblFATAL; } 
  /// the string name of exception category
  static string name(const log_exception &signal){ return "undefined exception";}
};

/// log_exceptions act as themselves
template<>
struct log_exception_traits<log_exception>{
  static int level(const log_exception &signal){ return log_exception::level(signal); } 
  static string name(const log_exception &signal){ return log_exception::name(signal); } 
};

# endif


# endif
