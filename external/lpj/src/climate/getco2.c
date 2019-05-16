/***************************************************************************/
/**                                                                       **/
/**                      g  e  t  c  o  2  .  c                           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  26.05.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define co2_p 280.0 /* preindustrial value for CO2 (ppm) */

Real getco2(const Climate *climate,int year)
{
  year-=climate->firstyear;
  return (year<0) ? co2_p : climate->co2[year];
} /* of 'getco2' */
