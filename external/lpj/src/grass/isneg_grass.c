/***************************************************************************/
/**                                                                       **/
/**               i  s  n  e  g  _  g  r  a  s  s  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.04.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

Bool isneg_grass(const Pft *pft, Real bm_inc)
{
  Pftgrass *grass;
  grass=pft->data;
#ifdef DEBUG2
  printf("isneg: %g %g\n",grass->ind.leaf,grass->ind.root);
#endif
  return (grass->ind.leaf<0 || grass->ind.root<0.0 || pft->fpc<=1e-20);
} /* of 'isneg_grass' */
