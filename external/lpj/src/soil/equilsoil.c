/***************************************************************************/
/**                                                                       **/
/**                     e  q  u  i  l  s  o  i  l  . c                    **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 20.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

/*
 *  DESCRIPTION
 *
 *  Analytically solves differential flux equations for fast and slow SOM 
 *  pools assuming annual litter inputs close to long term equilibrium
 *
 */

void equilsoil(Soil *soil) 
{
  
/*  soil->decomp_litter_mean/=soil_equil_year;
  soil->k_mean.fast/=soil_equil_year;
  soil->k_mean.slow/=soil_equil_year;*/
  soil->cpool.fast=soilfrac*fastfrac*soil->decomp_litter_mean/
                   soil->k_mean.fast;
  soil->cpool.slow=soilfrac*(1.0-fastfrac)*soil->decomp_litter_mean/
                   soil->k_mean.slow;
} /* of 'equilsoil' */
