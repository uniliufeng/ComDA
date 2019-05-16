/***************************************************************************/
/**                                                                       **/
/**              w  a  t  e  r  b  a  l  a  n  c  e  .  c                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.03.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define SOILDEPTH_EVAP 200.0  /* depth of sublayer at top of upper soil layer,
                                 from which evaporation is possible*/
#define BASEFLOW_FRAC 0.5     /* Fraction of standard percolation amount from
                                 lower soil layer*/
#define K_DEPTH 0.4 
#define K_AET 0.52            /* Fraction of total (vegetation) AET from upper
                                 soil layer*/
#define SOILDEPTH_UPPER soildepth[0]
#define PRIESTLEY_TAYLOR 1.32
#define  K_AET_DEPTH ((SOILDEPTH_UPPER/SOILDEPTH_EVAP-1.0)*(K_AET/K_DEPTH-1.0)/(1.0/K_DEPTH-1.0)+1.0)

Real waterbalance(Soil *soil,                /* Pointer to soil */
                  Real aet_stand[NSOILLAYER],
                  Real *evap,
                  Real rainmelt,             /* rainfall + melting water (mm) */
                  Real pet,
                  Real cover
                 )                           /* returns water runoff */
{
  Real runoff;
  Real perc,perc_frac,perc_baseflow;
  Real aet_evap;
  Real influx;
  int l;


 /* Weighting coefficient for AET flux from evaporation layer, assuming active
  *   root density decreases with soil depth
  *   Equates to 1.3 given SOILDEPTH_EVAP=200 mm, SOILDEPTH_UPPER=500 mm,
  *   K_DEPTH=0.4, K_AET=0.52 
  */

  runoff=0;
  if (cover>=1.0) cover=0.999;
  aet_evap=(aet_stand[0]*SOILDEPTH_EVAP*K_AET_DEPTH/SOILDEPTH_UPPER);
  *evap=pet*PRIESTLEY_TAYLOR*soil->w_evap*(1-cover);
  soil->w[0]+=(rainmelt-aet_stand[0]-*evap)/soil->par->whcs[0];
  /*printf("evap=%g cover=%g pet=%g w=%g\n",*evap,cover,pet,soil->w_evap);*/
#ifdef SAFE
  if (soil->w[0]<0.0){
    printf("rainmelt= %3.2f aet= %3.2f evap=  %3.2f cover=  %3.2f \n",rainmelt,aet_stand[0],*evap,cover);
    fail("Soil-moisture 1 negative: %g\n",soil->w[0]);
  }
#endif

  /*printf("rainmelt= %3.2f aet= %3.2f evap= %3.2f cover= %3.2f aet_evap= %3.2f w_evap= %3.2f w1= %3.2f\n",
   rainmelt,aet_stand[0],*evap,cover,aet_evap,soil->w_evap,soil->w[0]);*/

  if (soil->w[0]>1.0) 
  {
    /* Surface runoff */
    runoff=(soil->w[0]-1.0)*soil->par->whcs[0];
    soil->w[0]=1.0;
  }
  
  influx=rainmelt;
  
   /* Update water content in evaporation layer for tomorrow */

  soil->w_evap+=((rainmelt-*evap-aet_evap)/(soil->par->whc[0]*SOILDEPTH_EVAP));
  
  if (soil->w_evap>1.0) soil->w_evap=1.0;
  

   /* Percolation from evaporation layer */

  if (influx>0.1)
  {
    perc=SOILDEPTH_EVAP/soildepth[0]*soil->par->k1*
         pow(soil->w_evap,soil->par->k2);
    if(perc>influx)
      perc=influx;
    soil->w_evap-=perc/(soil->par->whc[0]*SOILDEPTH_EVAP);
  }
  
  if (soil->w_evap*(soil->par->whc[0]*SOILDEPTH_EVAP)>soil->w[0]*soil->par->whcs[0]) 
    soil->w_evap=soil->w[0]*soildepth[0]/SOILDEPTH_EVAP;
  if (soil->w_evap<0)
    soil->w_evap=0;


   /* Percolation and fluxes to and from lower soil layer(s)*/


  for (l=1;l<NSOILLAYER;l++) 
  {

    soil->w[l]-=aet_stand[l]/soil->par->whcs[l];
#ifdef SAFE 
    if (soil->w[l]<0.0)	     
     fail("Soil-moisture negative: %5.2f Soillayer: %4.3f\n",soil->w[l],l);
#endif
   
    /* Percolation: Allow only on days with rain or snowmelt */

    if (influx>=0.1)
    {
      perc=soil->par->k1*pow(soil->w[l-1],soil->par->k2);
      if(perc>influx)
        perc=influx;
      perc_frac=perc/soil->par->whcs[l-1];
      if(perc_frac>soil->w[l-1])
        perc_frac=soil->w[l-1];
      soil->w[l-1]-=perc_frac;
      soil->w[l]+=perc_frac*soil->par->whcs[l-1]/
                  soil->par->whcs[l];
      influx-=perc;
    }
 
    if (soil->w[l]>1.0) 
    {
      runoff+=(soil->w[l]-1.0)*soil->par->whcs[l];
      soil->w[l]=1.0;
    }

  }   
 
  /* Baseflow runoff (rain or snowmelt days only) */

  if (influx>=0.1) 
  {
    perc_baseflow=BASEFLOW_FRAC*soil->par->k1*pow(soil->w[BOTTOMLAYER],
                  soil->par->k2);
    if (perc_baseflow>influx) 
      perc_baseflow=influx;
    perc_frac=perc_baseflow/soil->par->whcs[BOTTOMLAYER];
    if(perc_frac>soil->w[BOTTOMLAYER])
      perc_frac=soil->w[BOTTOMLAYER];
    soil->w[BOTTOMLAYER]-=perc_frac;
    runoff+=perc_frac*soil->par->whcs[BOTTOMLAYER];
  }
  soil->meanw1+=soil->w[0];
  return runoff;
  
  /*Kf werte = k1*soildepth, wieviele mm pro tag können durch layer*/

} /* of 'waterbalance' */
