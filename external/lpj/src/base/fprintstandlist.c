/***************************************************************************/
/**                                                                       **/
/**               f  p  r  i  n  t  s  t  a  n  d  l  i  s  t  .  c       **/
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

void fprintstandlist(FILE *file,const Standlist standlist)
{
  Stand *stand;
  int s;
  fprintf(file,"Number of stands: %d\n",standlist->n);
  foreachstand(stand,s,standlist)
  {
    fprintf(file,"Stand: %d\n",s);
    fprintstand(file,stand);
  }
} /* of 'fprintstandlist' */
