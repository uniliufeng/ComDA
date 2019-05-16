/***************************************************************************/
/**                                                                       **/
/**                       i  n  i  t  p  f  t  .  c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 10.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/
#include "lpj.h"

void initpft(Pft *pft)
{
  pft->bm_inc=pft->wscal_mean=pft->phen=0;
  pft->par->init(pft);
} /* of 'initpft' */
