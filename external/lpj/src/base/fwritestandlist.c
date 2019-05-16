/***************************************************************************/
/**                                                                       **/
/**               f  w  r  i  t  e  s  t  a  n  d  l  i  s  t  .  c       **/
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

int fwritestandlist(FILE *file,const Standlist standlist,Bool full)
{
  Stand *stand;
  int s;
  fwrite(&standlist->n,sizeof(int),1,file);
  foreachstand(stand,s,standlist)
    if(fwritestand(file,stand,full))
      return s;
  return s;
} /* of 'fwritestandlist' */
