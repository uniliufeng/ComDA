/***************************************************************************/
/**                                                                       **/
/**              i  n  i  t  s  o  i  l  .  c                             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void initsoil(Soil *soil,const Soilpar *soilpar)
{
  int l;
  soil->par=soilpar;
  soil->cpool.fast=soil->cpool.slow=soil->litter.ag_tree=soil->litter.bg=
                   soil->litter.ag_grass=soil->k_mean.fast=soil->k_mean.slow=0.0;
  soil->alag=soil->amp=soil->meanw1=soil->decomp_litter_mean=0.0;
  for (l=0;l<NSOILLAYER;l++) 
    soil->w[l]=1.0;
  soil->w_evap=1.0;
} /* of 'initsoil' */
