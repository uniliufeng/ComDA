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
/**     Last change: 31.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void update_daily(Cell *cell, /* cell                   */
                  Real co2,   /* atmospheric CO2 (ppmv) */
                  Real temp,  /* temperature (deg C)    */
                  Real prec,  /* precipitation (mm)     */
                  Real sun,   /* sunshine (%)           */
                  int day    /* day (1..365)           */
                 )
{
  int p,s,l;
  Pft *pft;
  Real melt,pet,par,daylength,assim,evap,wet,rd,gpp,phen;
  Real gp_stand,gp_stand_leafon,cover_stand,cover_old;
  Real intercep_stand,fpc_total_stand;
  Real aet_stand[NSOILLAYER];
  Real gtemp_air,gtemp_soil;
  Stand *stand;
  petpar(&daylength,&par,&pet,cell->coord.lat,day,temp,sun);
#ifdef DEBUG3
      printf( "par= %.2f  pet= %.5f daylength= %.5f \n",par,pet,daylength);
      printf( "temp= %.2f sun= %.2f \n",temp,sun);
#endif

  gtemp_air=temp_response(temp);
  daily_climbuf(&cell->climbuf,temp);
  foreachstand(stand,s,cell->standlist)
  {
    /* allocate temporary arrays for phen and wet */
    foreachsoillayer(l)  aet_stand[l]=0.0;
    cover_stand=intercep_stand=0.0;
    gtemp_soil=temp_response(soiltemp(&stand->soil,&cell->climbuf));
    cell->output.mrh+=littersom(&stand->soil,gtemp_soil); 
    snow(&stand->soil.snowpack,&prec,&melt,temp);
    if(iffire && temp>0)
       stand->fire_sum+=fire_sum(&stand->soil.litter,stand->soil.w[0]);
    gp_stand=gp_sum(stand->pftlist,co2,temp,par,daylength,&gp_stand_leafon,&fpc_total_stand);
    foreachpft(pft,p,stand->pftlist)
    {
      phen=leaf_phenology(pft,temp,cell->coord.lat,day);
      cover_stand+=pft->fpc*phen;
      
      
      intercep_stand+=interception(&wet,pft,pet,prec);
/*      gp_stand=gp(pft,co2,temp,par,daylength,&gp_stand_leafon,&fpc_total_stand);*/
/*
 *    Calculate net assimilation, i.e. gross primary production minus leaf
 *    respiration, including conversion from FPC to grid cell basis.
 *
 */ 
      gpp=water_stressed(pft,aet_stand,stand->soil.w,gp_stand,gp_stand_leafon,fpc_total_stand,&rd,
                         &wet,pet,co2,temp,par,daylength);
      cell->output.mnpp+=npp(pft,pft->phen,gtemp_air,
                                          gtemp_soil,(gpp-rd));
                                              
#ifdef DEBUG3
      printf("PFT:%s fpc=%g aet1=%g aet2=%g phen=%g\n",pft->par->name,pft->fpc,aet_stand[0],aet_stand[1],phen);
      printf("PFT:%s gp_stand=%g phen=%g wscal=%g\n",pft->par->name,gp_stand,phen,pft->wscal);
#endif
                                              
    } /* of foreachpft */
    
    cell->output.mrunoff+=waterbalance(&stand->soil,aet_stand,&evap,
                                              prec+melt-intercep_stand,pet,
                                              cover_stand);
    cell->output.mevap+=evap; 
    foreachsoillayer(l)
    {
      cell->output.mtransp+=aet_stand[l];
      cell->output.mswc[l]+=stand->soil.w[l];
    }
    cell->output.minterc+=intercep_stand;
  } /* of foreachstand */
} /* of 'update_daily' */
