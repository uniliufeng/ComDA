/***************************************************************************/
/**                                                                       **/
/**              f  r  e  a  d  s  t  a  n  d  .  c                       **/
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

Stand *freadstand(FILE *file,const Pftpar pftpar[],const Soilpar *soilpar,
                  Bool swap)
{
  Stand *stand;
  stand=new(Stand);
  stand->frac=1.0;
  stand->pftlist=freadpftlist(file,pftpar,swap);
  freadsoil(file,&stand->soil,soilpar,swap);
  return stand;
} /* of 'freadstand' */
