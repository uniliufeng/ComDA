/***************************************************************************/
/**                                                                       **/
/**             w  a  t  e  r  _  s  t  r  e  s  s  e  d  .  c            **/
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

#define GM 3.26      /* empirical parameter in demand function */
#define EPSILON 0.1  /* min precision of solution in bisection method */
#define LAMBDAM 0.8  /* optimal Ci/Ca ratio */
#define ALPHAM 1.391 /* Priestley-Taylor coefficient (Demand,maximum-value)*/ 

typedef struct
{
  Real fac,co2,temp,par,daylength,tstress,fpar;
  int path;
} Data;

static Real fcn(Real lambda,Data *data)
{
  Real agd,rd;

/*
 *              Call photosynthesis to determine alternative total
 *              daytime photosynthesis estimate (adt2) implied by
 *              Eqns 2 & 19, Haxeltine & Prentice 1996, and current
 *              guess for lambda (xmid)
 */
  
  return data->fac*(1-lambda)-photosynthesis(&agd,&rd,data->path,
         lambda,data->tstress,data->co2,data->temp,data->par,data->fpar,
         data->daylength); 

/*
 *              Calculate total daytime photosynthesis implied by
 *              canopy conductance from water balance routine and
 *              current guess for lambda (xmid).  Units are mm/m2/day
 *              (mm come from gpd value, mm/day)
 *              Eqn 18, Haxeltine & Prentice 1996
 */

} /* of 'fcn' */

Real water_stressed(Pft *pft,
                    Real aet_layer[NSOILLAYER],
                    const Real w[NSOILLAYER], /* water in soil layers */
                    Real gp_stand,
                    Real gp_stand_leafon,
                    Real fpc,
                    Real *rd,
                    Real *wet,
                    Real pet,
                    Real co2,  /* Atmospheric CO2 partial pressure (ppm) */
                    Real temp, /* Temperature (deg C) */
                    Real par,
                    Real daylength /* Daylength (h) */
                   ) 
{
  int l; 
  Real supply, demand,wr,lambda,gpd,agd,gc,aet;
  Data data;
   
  wr=gpd=0;
  for (l=0;l<NSOILLAYER;l++) 
    wr+=pft->par->rootdist[l]*w[l];

  
  supply=pft->par->emax*wr*pft->phen;  		  /*supply=pft->par->emax*wr;*/
  /*demand=(pft->phen>0 && pft->fpc>0) ? (1.0-*wet)*pet*ALPHAM/(1+GM/((gp_stand+pft->par->gmin*pft->fpc)*pft->phen)) : 0;*/
  demand=(gp_stand>0) ? (1.0-*wet)*pet*ALPHAM/(1+GM/gp_stand) : 0;
  if  ( pet>0 && gp_stand_leafon>0 && pft->fpc>0)
  {
    /*pft->wscal=(pft->par->emax*wr*pft->fpc)/(pet*ALPHAM/(1+GM/(gp_stand_leafon+pft->par->gmin*pft->fpc)));*/
    pft->wscal=(pft->par->emax*wr)/(pet*ALPHAM/(1+GM/gp_stand_leafon));
    if(pft->wscal>1)
      pft->wscal=1;
  }
  else 
    pft->wscal=1;
  pft->wscal_mean+=pft->wscal;
  
  if(supply>=demand) 
    gc=gp_stand;
  else if(pet>0)
  {
    gc=GM*supply/((1.0-*wet)*pet*ALPHAM-supply);
    if(gc<0)
      gc=0;
  }
  else
    gc=0;
 
  aet=min(supply,demand)/wr*pft->fpc;
  for (l=0;l<NSOILLAYER;l++)
    aet_layer[l]+=aet*pft->par->rootdist[l]*w[l];
  
  /*gpd=hour2sec(daylength)*(gc-pft->par->gmin*pft->fpc*pft->phen);*/
  gpd=hour2sec(daylength)*(gc-pft->par->gmin*pft->phen)*pft->fpc;
  data.tstress=temp_stress(pft->par,temp,daylength);
  if(gpd>1e-5 && isphoto(data.tstress))
  {
    data.fac=gpd/1.6*ppm2bar(co2);
    data.path=pft->par->path;
    data.temp=temp;
    data.co2=ppm2Pa(co2);
    data.par=par;
    data.daylength=daylength;
    data.fpar=pft->phen*pft->fpc;
    lambda=bisect((Bisectfcn)fcn,0.02,LAMBDAM+0.05,&data,0,EPSILON,10); 
    photosynthesis(&agd,rd,data.path,lambda,data.tstress,data.co2,
                   temp,par,data.fpar,daylength);
    *rd=*rd;               
  }
  else
  {
    *rd=0; 
    return 0.0;
  }
  return agd;
} /* of 'water_stressed' */
