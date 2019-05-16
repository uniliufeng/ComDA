/***************************************************************************/
/**                                                                       **/
/**                u  p  d  a  t  e  _  d  a  i  l  y  .  c               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 18.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void update_daily(Cell *cell, /* cell                   */
                  Real co2,   /* atmospheric CO2 (ppmv) */
                  Real temp,  /* temperature (deg C)    */
                  Real prec,  /* precipitation (mm)     */
                  Real sun,   /* sunshine (%)           */
                  int day,    /* day (1..365)           */
                  int month   /* month (0..11)          */
                 )
{
  int p,s,l;
  Pft *pft;
  Real melt,pet,par,daylength,assim,evap,*phen,wet;
  Real gp_stand,cover_stand,intercep_stand;
  Real aet_stand[NSOILLAYER];
  Real gtemp_air,gtemp_soil;
  Stand *stand;
  petpar(&daylength,&par,&pet,cell->coord.lat,day,temp,sun);
  gtemp_air=temp_response(temp);
  daily_climbuf(&cell->climbuf,temp); 
  foreachstand(stand,s,cell->standlist)
  {
    /* allocate temporary arrays for phen and wet */
    phen=newvec(Real,getnpft(stand->pftlist));
    for (l=0;l<NSOILLAYER;l++) 
      aet_stand[l]=0.0;
    cover_stand=0.0;
    gtemp_soil=temp_response(soiltemp(&stand->soil,&cell->climbuf));
    cell->output.mrh[month]+=littersom(&stand->soil,gtemp_soil); 
    snow(&stand->soil.snowpack,&prec,&melt,temp);
    if(iffire && temp>0)
       stand->fire_sum+=fire_sum(&stand->soil.litter,stand->soil.w[0]);
    foreachpft(pft,p,stand->pftlist)
    {
      /* update generic part */
      phen[p]=leaf_phenology(pft,temp,cell->climbuf.gdd5,
                                    cell->coord.lat,day);
      cover_stand+=pft->fpc*phen[p];
    }
    gp_stand=gp_sum(stand->pftlist,co2,temp,par,daylength);
    intercep_stand=0;
    foreachpft(pft,p,stand->pftlist)
    {
/*
 *    Calculate net assimilation, i.e. gross primary production minus leaf
 *    respiration, including conversion from FPC to grid cell basis.
 *
 */ 
      intercep_stand+=interception(&wet,phen[p],pft,pet,prec);
      assim=water_stressed(pft,aet_stand,stand->soil.w,gp_stand,
                           phen[p],wet,pet,co2,temp,par,
                           daylength);
#ifdef DEBUG2
      printf( "pftvar[p].phen= %.2f \n assim= %.5f \n",pftvar[p].phen,assim);
      printf( "pftvar[p].wet= %.2f \n pet= %.2f \n",pftvar[p].wet,pet);
#endif
      cell->output.mnpp[month]+=npp(pft,phen[p],gtemp_air,
                                              gtemp_soil,assim);
    } /* of foreachpft */
    cell->output.mrunoff[month]+=waterbalance(&stand->soil,aet_stand,&evap,
                                              prec+melt-intercep_stand,pet,
                                              cover_stand);   
    cell->output.mevap[month]+=evap; 
    for(l=0;l<NSOILLAYER;l++)
    {
      cell->output.mtransp[month]+=aet_stand[l];
      cell->output.mswc[l][month]+=stand->soil.w[l];
    }
    cell->output.minterc[month]+=intercep_stand;
    free(phen);
  } /* of foreachstand */
} /* of 'update_daily' */
