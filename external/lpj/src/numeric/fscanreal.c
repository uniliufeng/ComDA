/***************************************************************************/
/**                                                                       **/
/**                   f  s  c  a  n  r  e  a  l  .  c                     **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"

Bool fscanreal(FILE *file,Real *val)
{
  double x;
  if(fscanf(file,"%lg",&x)!=1)
    return TRUE;
  *val=x;
  return FALSE;
} /* of 'fscanreal' */
