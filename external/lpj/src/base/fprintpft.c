/***************************************************************************/
/**                                                                       **/
/**                 f  p  r  i  n  t  p  f  t  .  c                       **/
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


void fprintpft(FILE *file,const Pft *pft)
{
  fprintf(file,"Pft: %s\n",pft->par->name);
  fprintf(file,"fpc: %g\n",pft->fpc);
  fprintf(file,"nind: %g\n",pft->nind);
  fprintf(file,"wscal: %g\n",pft->wscal_mean);
  fprintf(file,"aphen: %g\n",pft->aphen);
  fprintf(file,"bminc: %g\n",pft->bm_inc);
  fprintf(file,"gdd: %g\n",pft->gdd);
  pft->par->fprint(file,pft);
} /* of 'fwritepft' */
