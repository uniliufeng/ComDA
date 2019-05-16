/***************************************************************************/
/**                                                                       **/
/**                      f  a  i  l  .  c                                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "types.h"

void fail(const char *msg,...)
{
  va_list ap;
  va_start(ap,msg);
  vfprintf(stderr,msg,ap);
  fprintf(stderr,", program terminated unsuccessfully.\n");
  va_end(ap);
  exit(EXIT_FAILURE);
} /* of 'fail' */
