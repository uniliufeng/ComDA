/***************************************************************************/
/**                                                                       **/
/**                 f  p  r  i  n  t  s  t  a  n  d  .  c                 **/
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

#include "lpj.h"

void fprintstand(FILE *file,const Stand *stand)
{
  fprintsoil(file,&stand->soil);
  fprintpftlist(file,stand->pftlist);
} /* of 'fprintstand' */
