/***************************************************************************/
/**                                                                       **/
/**                       n  e  w  p  f  t  .  c                          **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 24.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Pft *newpft(const Pftpar *pftpar /* Parameter of pft */
           )                     /* returns allocated pointer to pft */
{
  Pft *pft;
  pft=new(Pft);
  pft->par=pftpar;
  pft->gddbuf=newbuffer(CLIMBUFSIZE);
  pft->fpc=pft->nind=pft->wscal=pft->aphen=pft->bm_inc=pft->wscal_mean=pft->gdd=0.0;
  pft->par->newpft(pft);
  return pft;
} /* of 'newpft' */
