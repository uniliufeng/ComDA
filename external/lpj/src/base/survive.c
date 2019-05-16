/***************************************************************************/
/**                                                                       **/
/**                  s  u  r  v  i  v  e  .  c                            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 04.05.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool survive(const Pftpar *pftpar,const Climbuf *climbuf)
{
  Real temp_min20,temp_max20;
  
  temp_min20=climbuf->min->avg/climbuf->min->n;
  temp_max20=climbuf->max->avg/climbuf->max->n;
  return (temp_min20>=pftpar->temp.low) || 
         (temp_max20-temp_min20>=pftpar->min_temprange);

} /* of 'survive' */
