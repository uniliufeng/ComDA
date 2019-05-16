/***************************************************************************/
/**                                                                       **/
/**                      l  i  t  t  e  r  s  o  m  .  c                  **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 28.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define MOIST_DENOM (1.0-exp(-1.0))
#define k_litter10 (0.3/NDAYYEAR)
#define k_soilfast10 (0.03/NDAYYEAR)
#define k_soilslow10 (0.001/NDAYYEAR)
#define cpool_low 1e-5

Real littersom(Soil *soil,Real gtemp_soil)
{
  Real response,decay_litter,cflux_litter,cflux_soil_fast,cflux_soil_slow;
  Litter decom;
  if(gtemp_soil>0)
  {
    response=gtemp_soil*((1.0-exp(-soil->w[0]))/MOIST_DENOM);

/*
 *       Calculate monthly decomposition rates (k, /month) as a function of
 *       temperature and moisture
 *
 */

    decay_litter=1.0-exp(-(k_litter10*response));
   /* decay_litter=k_litter10*response;*/
    decom.ag_tree=soil->litter.ag_tree*decay_litter;
    decom.ag_grass=soil->litter.ag_grass*decay_litter;
    decom.bg=soil->litter.bg*decay_litter;
    soil->litter.ag_tree-=decom.ag_tree;
    soil->litter.ag_grass-=decom.ag_grass;
    soil->litter.bg-=decom.bg;
    cflux_litter=soilfrac*littersum(decom);
    soil->cpool.fast+=fastfrac*cflux_litter;
    soil->cpool.slow+=(1-fastfrac)*cflux_litter;
    cflux_soil_slow=soil->cpool.slow*(1.0-exp(-(k_soilslow10*response)));
    cflux_soil_fast=soil->cpool.fast*(1.0-exp(-(k_soilfast10*response)));
    /*cflux_soil_slow=soil->cpool.slow*k_soilslow10*response;
    cflux_soil_fast=soil->cpool.fast*k_soilfast10*response;*/
    soil->cpool.slow-=cflux_soil_slow;
    soil->cpool.fast-=cflux_soil_fast; 
    /*sum for equilsom-routine*/
    soil->decomp_litter_mean+=littersum(decom)/soil_equil_year;
    soil->k_mean.fast+=(k_soilfast10*response)/soil_equil_year;
    soil->k_mean.slow+=(k_soilslow10*response)/soil_equil_year;
    /* Empty soil pools below a minimum threshold */
    setminthreshold(soil->cpool.fast,cpool_low);
    setminthreshold(soil->cpool.slow,cpool_low);
    return cflux_litter*atmfrac/soilfrac+cflux_soil_slow+cflux_soil_fast;
  }
  else
    return 0.0;
} /* of 'littersom' */
