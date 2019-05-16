/***************************************************************************/
/**                                                                       **/
/**               f  r  e  a  d  p  f  t  .  c                            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 24.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Pft *freadpft(FILE *file,const Pftpar pftpar[],Bool swap)
{
  Pft *pft;
  int id;
  if(freadint1(&id,swap,file)!=1)
    return NULL;
  pft=new(Pft);
  pft->par=pftpar+id;
  /*fread(&pft->fpc,1,swap,file);*/
  /*fread(&pft->fpar,,1,swap,file);*/
  freadreal1(&pft->wscal,swap,file);
  freadreal1(&pft->wscal_mean,swap,file);
  freadreal1(&pft->aphen,swap,file);
  freadreal1(&pft->phen,swap,file);
  pft->par->fread(file,pft,swap);
  freadreal1(&pft->bm_inc,swap,file);
  freadreal1(&pft->nind,swap,file);
  freadreal1(&pft->gdd,swap,file);
  freadreal1(&pft->fpc,swap,file);
  pft->gddbuf=freadbuffer(file,swap);
  return pft;
} /* of 'freadpft' */
