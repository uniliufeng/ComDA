/***************************************************************************/
/**                                                                       **/
/**             p  h  e  n  o  l  o  g  y  _  g  r  a  s  s  .  c         **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

#define COLDEST_DAY_NHEMISPHERE 14
#define COLDEST_DAY_SHEMISPHERE 195

Real phenology_grass(Pft *pft,Real temp,Real lat,int day)
{
  Real phen,dtemp;
  dtemp=temp - getpftpar(pft,gddbase);
  if(dtemp>0.0)
    pft->gdd+=dtemp;
  pft->phen=pft->gdd*pft->par->ramp;
  if(pft->phen>1)
    pft->phen=1;
  if ((lat>=0.0 && day==COLDEST_DAY_NHEMISPHERE) ||
      (lat<0.0 && day==COLDEST_DAY_SHEMISPHERE)) 
    pft->aphen=pft->gdd=pft->phen=0.0;
  pft->aphen+=pft->phen;
  phen=pft->phen;
  return phen;
} /* of 'phenology_grass' */
