/*e***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD
 *
 *   $Revision: 1.16 $
 *   $Date: 2013/08/16 11:32:45 $
 *   @(#) $Header: /home/plasmacvs/source_tree/gridmd/include/wxinc.h,v 1.16 2013/08/16 11:32:45 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/gridmd/include/wxinc.h,v $
$Revision: 1.16 $
$Author: valuev $
$Date: 2013/08/16 11:32:45 $
*/
/*e****************************************************************************
 * $Log: wxinc.h,v $
 * Revision 1.16  2013/08/16 11:32:45  valuev
 * compiled tcpengine
 *
 * Revision 1.15  2013/03/07 17:08:59  valuev
 * kintech xml
 *
 * Revision 1.14  2012/11/22 17:01:06  valuev
 * sync with kintech svn
 *
 * Revision 1.13  2012/09/07 14:26:06  valuev
 * removed using namespace from h files
 *
 * Revision 1.12  2012/09/04 13:14:05  valuev
 * wx compatibility
 *
 * Revision 1.11  2012/08/31 15:15:44  valuev
 * restructured gridmd workflow
 *
 * Revision 1.10  2011/07/22 14:24:55  valuev
 * forces with norm matrix
 *
 * Revision 1.9  2011/07/21 14:33:16  valuev
 * made compatible with wxWidgets 2.9 (unicode)
 *
 * Revision 1.8  2010/10/07 11:20:31  valuev
 * preliminary program restart
 *
 * Revision 1.7  2010/10/01 07:01:04  valuev
 * added xml resource saving
 *
 * Revision 1.6  2010/06/12 17:57:31  valuev
 * some workflow coding
 *
 * Revision 1.5  2009/06/10 20:53:28  valuev
 * updated splindes, trajReader
 *
 * Revision 1.4  2008/02/21 14:02:41  valuev
 * Added parametric methods
 *
 * Revision 1.3  2007/12/13 18:58:38  valuev
 * Added interface to Dakota
 *
 * Revision 1.2  2007/11/20 22:15:48  valuev
 * Added parametric project
 *
 * Revision 1.1  2005/12/02 18:50:56  valuev
 * added  HEAD project tree
 *
 * Revision 1.1  2005/11/30 23:37:39  valuev
 * put gridmd on cvs at biolab1.mipt.ru
 *
 *
*******************************************************************************/
# ifndef __WXINC_H
# define __WXINC_H

/*e @file wxinc.h
    @brief Necessary includes/definitions for wxWidgets Base library    
**/

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "wx/defs.h"

#if wxUSE_GUI
    #error "This project can't be compiled in GUI mode."
#endif // wxUSE_GUI
#include <stdio.h>

//e default is STL_HASH
# ifndef STL_HASH
# define STL_HASH 1
# endif

//e These macros are used to put wxWidgets declarations into separate namespace.
//e Trick for VC class browser only.



#include "wx/string.h"
#include "wx/file.h"
#include "wx/log.h"
#include "wx/app.h"

#include "wx/filename.h"
#include "wx/tokenzr.h"
//#include "wx/utils.h"
//#include "wx/process.h"
# include "wx/cmdline.h"

#ifndef NO_XML
#include "wx/xml/xml.h"
#endif

# ifndef STL_HASH
# define gm_DECLARE_STRING_HASH_MAP(type,hash_type)  \
         namespace wxdef { \
         WX_DECLARE_STRING_HASH_MAP(type,hash_type); \
         };\
         using namespace wxdef;
# else
# include <map>
# define gm_DECLARE_STRING_HASH_MAP(type,hash_type) typedef std::map<wxString,type> hash_type;
# endif



    // Single inheritance with one base class
#define IMPLEMENT_DYNAMIC_CLASS_12(name, basename, firstname)                               \
    wxIMPLEMENT_CLASS_COMMON1(name, basename, name::wxCreateObject)           \
    wxObject* name::wxCreateObject()                                          \
        { return (firstname *)new name; }




// without this pragma, the stupid compiler precompiles #defines below so that
// changing them doesn't "take place" later!
#ifdef __VISUALC__
    #pragma hdrstop
#endif


#if wxMAJOR_VERSION >= 3 

inline const char *fmt(const wxCStrData format,...){
  va_list args;
  va_start(args,format);
  static char buff[1024];
  vsnprintf(buff,1024,(const char *)format,args);
  return buff;
}


# define GetPropVal GetAttribute
# define DeleteProperty DeleteAttribute
# define AddProperty AddAttribute

# endif

# endif
