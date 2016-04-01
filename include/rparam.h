/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.7 $
 *   $Date: 2015/01/30 09:08:04 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/rparam.h,v 1.7 2015/01/30 09:08:04 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/rparam.h,v $
$Revision: 1.7 $
$Author: valuev $
$Date: 2015/01/30 09:08:04 $
*/
/*s****************************************************************************
 * $Log: rparam.h,v $
 * Revision 1.7  2015/01/30 09:08:04  valuev
 * working on sweeps
 *
 * Revision 1.6  2009/02/10 14:20:45  valuev
 * sync with FDTD project
 *
 * Revision 1.5  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.5  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.4  2007/02/20 10:26:11  valuev
 * added newlines at end of file
 *
 * Revision 1.3  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
 * Revision 1.2  2006/08/08 13:02:04  valuev
 * Added geometry
 *
 * Revision 1.2  2006/07/21 16:22:03  valuev
 * Added Tight Binding for graphite+O
 *
 * Revision 1.1  2005/12/02 18:51:06  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:36:11  valuev
 * put ivutils to cvs on biolab1.mipt.ru
 *
 * Revision 1.1  2005/11/30 23:15:43  valuev
 * put ivutils on cvs biolab1.mipt.ru
 *
 *
*******************************************************************************/
# ifndef __rparam_h
# define __rparam_h

# include <stdarg.h>
# include <stdio.h>

# ifdef __cplusplus
extern "C" {
# endif

int   Open_param_file(const char *file); // returns !0 if OK
int   Read_param(const char *format,...);
int   Read_paramn(int n,const char *format,...);
void  Close_param_file(void);
long  Set_position(const char *format);
void  Set_cur_pos(long lpos);
long  Get_cur_pos(void);
void  Set_comment_char(char c);
void  Set_help_char(char c);
char  Get_help_char();
char  Get_comment_char();
void  Set_wildcard_char(char c);
char  Get_wildcard_char();

                        // In case of errors:
void  Set_stop(int s);  // 1-all errors fatal,
                        // -1 -shows error messages and returns
                        // 0 - doesn't show any messages and returns
int   Get_stop(void);
char  Set_delimiter(char c);
char  *Get_delimiter();
FILE  *Get_cur_file(void);
char *Get_cur_filename(void);

void  Set_scan_end(char c);

int Get_cur_help(char *hbuf,int maxlen);

extern char AcFormat[];

# ifdef __cplusplus
}
# endif



# if 0 // moved to common.h
//# ifdef FUCK

int vsscanf(char *,char *, va_list);

# endif

/* extern int vsscanf(char * str,char *format,va_list args);*/



# endif













