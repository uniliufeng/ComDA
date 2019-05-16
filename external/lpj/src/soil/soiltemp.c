/***************************************************************************/
/**                                                                       **/
/**                    s  o  i  l  t  e  m  p  .  c                       **/
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

#define LAG_CONV (NDAYYEAR*0.5*M_1_PI)  /* conversion factor for oscillation 
                                lag from angular units to days (=365/(2*PI))*/

Real soiltemp(const Soil *soil,const Climbuf *climbuf)
{
  int i;
  Real a,b,temp_lag;
  if(soil->w[0]==0)
    return climbuf->temp[NDAYS-1];
  linreg(&a,&b,climbuf->temp,NDAYS);
  temp_lag=a+b*(NDAYS-1-soil->alag*LAG_CONV);
  return climbuf->atemp_mean+soil->amp*(temp_lag-climbuf->atemp_mean);
} /* of 'soiltemp' */
