/***************************************************************************/
/**                                                                       **/
/**           v  e  g  c  _  s  u  m  _  g  r  a  s  s  .  c              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

Real vegc_sum_grass(const Pft *pft)
{
  const Pftgrass *grass;
  grass=pft->data;
  return phys_sum_grass(grass->ind)*pft->nind;
} /* of 'vegc_sum_grass' */
