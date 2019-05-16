/***************************************************************************/
/**                                                                       **/
/**                f  i  r  e  _  s  u  m  .  c                           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 18.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define flam_tree 0.3
#define flam_grass 0.3

Real fire_sum(const Litter *litter,Real w_surf)
{
  Real moistfactor,litter_sum;
  litter_sum=litter->ag_tree+litter->ag_grass;
  if(litter_sum==0)
    return 0;
  moistfactor=(flam_tree*litter->ag_tree+flam_grass*litter->ag_grass)/
              litter_sum;
  return (moistfactor>0) ? 
            exp(-M_PI*(w_surf/moistfactor)*(w_surf/moistfactor)) : 0;
} /* of 'fire_sum' */
