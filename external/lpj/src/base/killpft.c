/***************************************************************************/
/**                                                                       **/
/**                         k  i  l  l  p  f  t  .  c                     **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  15.11.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool killpft(Litter *litter,Pft *pft,const Climbuf *climbuf)
{
  if(!survive(pft->par,climbuf))
  {
    litter_update(litter,pft,pft->nind);
    return TRUE; 
  }
  else
    return FALSE;
} /* of 'killpft' */
