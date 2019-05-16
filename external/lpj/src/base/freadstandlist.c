/***************************************************************************/
/**                                                                       **/
/**             f  r  e  a  d  s  t  a  n  d  l  i  s  t  .  c            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Standlist freadstandlist(FILE *file,const Pftpar pftpar[],
                         const Soilpar *soilpar,Bool swap)
{
  int i,n;
  Standlist standlist;
  if(freadint1(&n,swap,file)!=1)
    return NULL;
  standlist=newlist();
  for(i=0;i<n;i++)
    addlistitem(standlist,freadstand(file,pftpar,soilpar,swap));
  return standlist;
} /* of 'freadstandlist' */
