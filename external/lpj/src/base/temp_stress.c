/***************************************************************************/
/**                                                                       **/
/**                t  e  m  p  _  s  t  r  e  s  s  .  c                  **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 11.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define tmc3 45.0   /* maximum temperature for C3 photosynthesis */
#define tmc4 55.0   /* maximum temperature for C4 photosynthesis */

/*
 *     Function temp_stress for photosynthesis
 *
 *     Adapted from Farquhar (1982) photosynthesis model, as simplified by
 *     Collatz et al 1991, Collatz et al 1992 and Haxeltine & Prentice 1996
 *
 */

Real temp_stress(const Pftpar *pftpar,Real temp,Real daylength)
{
  Real k1,k2,k3,low,high;
  if(daylength<0.01 || (pftpar->path==C3 && temp>tmc3) 
                    || (pftpar->path==C4 && temp>tmc4))
    return 0.0;
  if(temp<pftpar->temp_co2.high)
  {
    k1=2*log(1/0.99-1)/(pftpar->temp_co2.low-pftpar->temp_photos.low);
    k2=(pftpar->temp_co2.low+pftpar->temp_photos.low)*0.5;
    low=1/(1+exp(k1*(k2-temp)));
    k3=log(0.99/0.01)/(pftpar->temp_co2.high-pftpar->temp_photos.high);
    high=1-0.01*exp(k3*(temp-pftpar->temp_photos.high));
    return low*high;
  }
  return 0.0;
} /* of 'temp_stress' */
