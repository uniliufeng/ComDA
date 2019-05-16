/***************************************************************************/
/**                                                                       **/
/**                f  i  r  e  _  p  r  o  b  .  c                        **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 09.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define minfuel 200.0 /* fuel threshold to carry fire (gC/m2) */

Real fire_prob(const Litter *litter,Real fire_sum,Real *fire_frac)
{
  Real fire_index,sm;
  fire_index=fire_sum/NDAYYEAR;
  sm=fire_index-1;
  *fire_frac=fire_index*exp(sm/(0.45*sm*sm*sm+2.83*sm*sm+2.96*sm+1.04));
#ifdef SAFE
  if (*fire_frac>1.0)
    fail("fire: probability of fire=%g >1.0",*fire_frac);
#endif
  return (*fire_frac<0.001 ||(litter->ag_tree+litter->ag_grass<minfuel)) ?
            0.001 : *fire_frac;
} /* of 'fire_prob' */
