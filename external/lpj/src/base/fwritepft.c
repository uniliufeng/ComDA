/***************************************************************************/
/**                                                                       **/
/**                 f  w  r  i  t  e  p  f  t  .  c                       **/
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


Bool fwritepft(FILE *file,const Pft *pft,Bool full)
{
  fwrite1(&pft->par->id,sizeof(int),file);
  /*fwrite1(&pft->fpc,sizeof(Real),file);*/
  /*ffwrite1(&pft->fpar,sizeof(Real),file);*/
  fwrite1(&pft->wscal,sizeof(Real),file);
  fwrite1(&pft->wscal_mean,sizeof(Real),file);
  fwrite1(&pft->aphen,sizeof(Real),file);
  fwrite1(&pft->phen,sizeof(Real),file);
  pft->par->fwrite(file,pft);
  if(full)
  {
    fwrite1(&pft->bm_inc,sizeof(Real),file);
    fwrite1(&pft->nind,sizeof(Real),file);
    fwrite1(&pft->gdd,sizeof(Real),file);
    fwrite1(&pft->fpc,sizeof(Real),file);
    if(fwritebuffer(file,pft->gddbuf))
      return TRUE;
  }
  return FALSE;
} /* of 'fwritepft' */
