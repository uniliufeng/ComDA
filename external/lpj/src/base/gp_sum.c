/***************************************************************************/
/**                                                                       **/
/**                     g  p  _  s  u  m  .  c                            **/
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

#define lambda_opt 0.8

Real gp_sum(const Pftlist pftlist, /* Pft list */
            Real co2,              /* atmospheric CO2 concentration (ppm) */
            Real temp,             /* temperature (deg C) */
            Real par,              /* photosynthetic active radiation flux */
            Real daylength,        /* daylength (h) */
            Real *gp_stand_leafon, /* pot. canopy conduct.at full leaf cover */
            Real *fpc		   /* total fpc of all Pfts */
           )
{
  int p;
  const Pft *pft;
  Real agd,adtmm,gp,rd,fpc_total;
  *gp_stand_leafon=gp=*fpc=0;
  if(daylength==0)
    return 0;
  foreachpft(pft,p,pftlist)
  {
    adtmm=photosynthesis(&agd,&rd,pft->par->path,lambda_opt,
                         temp_stress(pft->par,temp,daylength),ppm2Pa(co2),
                         temp,par,pft->fpc,daylength);
     gp+=((1.6*adtmm/(ppm2bar(co2)*(1.0-lambda_opt)*hour2sec(daylength)))+
         pft->par->gmin*pft->fpc)*pft->phen;
    *gp_stand_leafon+=(1.6*adtmm/(ppm2bar(co2)*(1.0-lambda_opt)*hour2sec(daylength)))+
                    pft->par->gmin*pft->fpc;
    *fpc+=pft->fpc;
  }
  fpc_total=*fpc;
  *gp_stand_leafon= (gp<1e-20 || fpc_total<1e-20) ? 0 : *gp_stand_leafon/fpc_total;
  return (gp<1e-20 || fpc_total<1e-20) ? 0 : gp/fpc_total;
} /* of 'gp_sum' */

