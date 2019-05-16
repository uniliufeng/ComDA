/***************************************************************************/
/**                                                                       **/
/**              e  s  t  a  b  l  i  s  h  .  c                          **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 03.05.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool establish(Real gdd,const Pftpar *pftpar,const Climbuf *climbuf)
{
  Real temp_min20;
  temp_min20=climbuf->min->avg/climbuf->min->n;
  return (temp_min20>=pftpar->temp.low) && 
         (temp_min20<=pftpar->temp.high) && 
         (gdd>=pftpar->gdd5min);
/*         (gdd>=pftpar->gdd5min) && (climbuf->temp_max<=pftpar->twmax);*/
} /* of 'establish' */
