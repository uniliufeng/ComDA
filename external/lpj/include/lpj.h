/***************************************************************************/
/**                                                                       **/
/**                      l  p  j  .  h                                    **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     LPJ header file contains all necessary includes                   **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef LPJ_H /* Already included? */
#define LPJ_H

#define LPJ_VERSION  "0.9.025"  /*  now beta, validation in progress */

/* Necessary header files */

/* Standard C header files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

/*  Defined header files for LPJ */

#include "list.h"
#include "types.h"
#include "popen.h"
#include "swap.h"
#include "coord.h"
#include "config.h"
#include "buffer.h"
#include "date.h"
#include "climate.h"
#include "soil.h"
#include "pft.h"
#include "errmsg.h"
#include "crop.h"
#include "pftlist.h"
#include "numeric.h"
#include "units.h"
#include "conf.h"
#include "stand.h"
#include "output.h"
#include "header.h"
#include "cell.h"

extern Bool iffire;

/* Declaration of functions */

extern Cell *newgrid(Config *,const Pftpar [],int,const Soilpar [],int);
extern Bool fwriterestart(const char *,const Cell *,int,int,int,int);

/* Definition of macros */

#define israndomprec(config) (config.wet_filename!=NULL)
#define iswriterestart(config) (config.write_restart_filename!=NULL)
#define isreadrestart(config) (config.restart_filename!=NULL)

#endif
