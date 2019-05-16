/***************************************************************************/
/**                                                                       **/
/**                    f  i  r  e  p  f  t  .  c                          **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Real firepft(Litter *litter,Pftlist pftlist,Real fire_frac)
{
  int p;
  Pft *pft;
  Real flux,flux_litter;
  flux=flux_litter=0;
  foreachpft(pft,p,pftlist)
    flux+=fire(pft,fire_frac);
    flux_litter=(litter->ag_grass+litter->ag_tree)*fire_frac;
  litter->ag_grass*=(1-fire_frac);
  litter->ag_tree*=(1-fire_frac);
  return flux+flux_litter; 
} /* of 'firepft' */
