/***************************************************************************/
/**                                                                       **/
/**                    g  e  t  l  a  g  .  c                             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 16.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define DEPTH (soildepth[0]*0.0005)  /*soil depth at which to estimate temperature (m)*/
#define DIFFUS_CONV 0.0864          /*Convert diffusivity from mm2/s to m2/day*/
#define HALF_OMEGA (M_PI/NDAYYEAR)     /*corresponds to omega/2 = pi/365*/

void getlag(Soil *soil,int month)
{
  Real diffus;
  soil->meanw1*=ndaymonth1[month];
  diffus=(soil->meanw1<0.15)  ?
            (soil->par->tdiff_15-soil->par->tdiff_0)/0.15*soil->meanw1+
            soil->par->tdiff_0 :
            (soil->par->tdiff_100-soil->par->tdiff_15)/0.85*(soil->meanw1-0.15)
            +soil->par->tdiff_15;
  soil->meanw1=0;
  soil->alag= DEPTH/sqrt(diffus*DIFFUS_CONV/HALF_OMEGA);
  soil->amp=exp(-soil->alag);
} /* of 'getlag' */
