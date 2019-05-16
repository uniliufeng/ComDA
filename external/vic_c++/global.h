/**********************************************************************
                        Global Variables

  NOTE: This file exists because global variables that are shared among
        files via the "extern" statement must be initially declared
        (without the word "extern") ONLY once.  Currently, vicNl_def.h
        is included (via vicNl.h) in every .c file, meaning that any
        declarations in vicNl_def.h end up happening multiple times
        (once per .c file).  Thus, these "extern" variables cannot be
        declared in vicNl_def.h.  This is not a problem for #define
        statements and typedef statements, which is what vicNl_def.h
        is primarily composed of.

  $Id: global.h,v 4.9.2.8 2009/10/20 18:42:07 vicadmin Exp $

  29-Oct-03 Added version string and removed unused options from
	    optstring.							TJB
  2009-Jun-09 Added definitions of reference landcover types, used
	      mainly for pot_evap computations but also defines the
	      characteristics of bare soil.				TJB
**********************************************************************/
#ifndef global_H_
#define global_H_
#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif  
#endif
