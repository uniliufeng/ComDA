/***************************************************************************/
/**                                                                       **/
/**                p  h  o  t  o  s  y  n  t  h  e  s  i  s  .  c         **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define po2 20.9e3/* O2 partial pressure in Pa */
#define p 1.0e5   /* atmospheric pressure in Pa */
#define bc3 0.015 /* leaf respiration as fraction of Vmax for C3 plants */
#define bc4 0.02  /* leaf respiration as fraction of vmax for C4 plants */
#define theta 0.7 /* colimitation (shape) parameter */
#define q10ko 1.2 /* q10 for temperature-sensitive parameter ko */
#define q10kc 2.1 /* q10 for temperature-sensitive parameter kc */
#define q10tau 0.57 	/* q10 for temperature-sensitive parameter tau */
#define ko25 3.0e4  	/* value of ko at 25 deg C */
#define kc25 30.0   	/* value of kc at 25 deg C */
#define tau25 2600.0	/* value of tau at 25 deg C */
#define alphaa 0.5  	/* fraction of PAR assimilated at ecosystem level */
                              /* relative to leaf level */
#define alphac3 0.08  /* intrinsic quantum efficiency of CO2 uptake in */
                                /* C3 plants */
#define alphac4 0.053 /* C4 intrinsic quantum efficiency */
#define cmass 12.0    /* atomic mass of carbon */
#define cq 4.6e-6   /* conversion factor for solar radiation at 550 nm */
                              /* from J/m2 to E/m2 (E mol quanta) */
#define lambdamc4 0.4/* optimal ratio of intercellular to ambient CO2 */
                               /* concentration (lambda) in C4 plants */
#define lambdamc3 0.8/* optimal (maximum) lambda in C3 plants */
#define n0 7.15     /* leaf N concentration (mg/g) not involved in */
                              /* photosynthesis */
#define m 25.0      /* corresponds to #define p in Eqn 28, Haxeltine */
                              /*  Prentice 1996 */
#define t0c3 250.0 /* base temperature (K) in Arrhenius temperature */
                              /* response function for C3 plants */
#define t0c4 260.0  /* base temperature in Arrhenius func for C4 plants */
#define e0 308.56   /* parameter in Arrhenius temp response function */
#define tk25 298.15 /* 25 deg C in Kelvin */

/*
 *     Function photosynthesis 
 *
 *     Adapted from Farquhar (1982) photosynthesis model, as simplified by
 *     Collatz et al 1991, Collatz et al 1992 and Haxeltine & Prentice 1996
 *
 */

Real photosynthesis(Real *agd,
		    Real *rd,
                    int path,       /* Path (C3/C4) */ 
                    Real lambda,
                    Real tstress,
                    Real co2,       /* atmospheric CO2 partial pressure (Pa) */
                    Real temp,      /* temperature (deg C) */
                    Real par,
                    Real fpar,
                    Real daylength  /* daylength (h) */
                   )
{
  static Real ko,kc,tau,pi,c1,c2;
  Real apar,je,jc,phipi,and,adt,b,s,sigma,vm;
  static Real fac,gammastar,temp_save=9999;
  if(tstress<1e-2) 
  {
    *agd=0;
    *rd=0;
    return 0;
  }
  else
  {
    apar=fpar*par*alphaa;
    if(path==C3)
    {
      if(temp!=temp_save)
      {
        temp_save=temp;
        ko=ko25*exp((log(q10ko))*(temp-25)*0.1);
        kc=kc25*exp((log(q10kc))*(temp-25)*0.1);
        fac=kc*(1+po2/ko);
        tau=tau25*exp((log(q10tau))*(temp-25)*0.1);
        gammastar=po2/(2*tau);
      }
      pi=lambdamc3*co2;
      c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar));
 
      /* Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996 */

      c2=(pi-gammastar)/(pi+fac);
      
      s=(24/daylength)*bc3;
      sigma=1-(c2-s)/(c2-theta*s);
      sigma= (sigma<=0) ? 0 : sqrt(sigma);
      b=bc3;   /* Choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996 */
      /*
       *       Intercellular CO2 partial pressure in Pa
       *       Eqn 7, Haxeltine & Prentice 1996
       */
      
      vm=(1.0/bc3)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*
         cmass*cq;

      pi=lambda*co2;

      /*       Recalculation of C1C3, C2C3 with actual pi */

      c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar));

      c2=(pi-gammastar)/(pi+fac);
    }
    else /* C4 photosynthesis */
    {
      temp_save=999;
      c1=tstress*alphac4;
      c2=1.0;
      b=bc4;
      s=(24/daylength)*bc4;
      sigma=1-(c2-s)/(c2-theta*s);
      sigma= (sigma<=0) ? 0 : sqrt(sigma);
      vm=(1.0/bc4)*c1/c2*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*
         cmass*cq;
      
      /*
       *       Parameter accounting for effect of reduced intercellular CO2
       *       concentration on photosynthesis, Phipi.
       *       Eqn 14,16, Haxeltine & Prentice 1996
       *       Fig 1b, Collatz et al 1992
       */
      phipi=lambda/lambdamc4;
      if(phipi<1)
        c1=tstress*phipi*alphac4;
    }
    
    /*
     *     je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
     *     Convert je from daytime to hourly basis
     *
     *     Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
     *     Eqn 3, Haxeltine & Prentice 1996
     */
    je=c1*apar*cmass*cq/daylength;

    /*
     * Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
     *    Eqn 5, Haxeltine & Prentice 1996
     */
    jc=c2*hour2day(vm);

    /*
     *    Calculation of daily gross photosynthesis, Agd, gC/m2/day
     *    Eqn 2, Haxeltine & Prentice 1996
     */
    
    *agd=(je+jc-sqrt((je+jc)*(je+jc)-4.0*theta*je*jc))/(2.0*theta)*daylength;

    /*    Daily leaf respiration, Rd, gC/m2/day
     *    Eqn 10, Haxeltine & Prentice 1996
     */
    *rd=b*vm;

    /*    Daily net photosynthesis (at leaf level), And, gC/m2/day */

    and=*agd-*rd;
  
    /*     Total daytime net photosynthesis, Adt, gC/m2/day
     *     Eqn 19, Haxeltine & Prentice 1996
     */
 
    adt=and+(1.0-hour2day(daylength))*(*rd); 
       
    /*     Convert adt from gC/m2/day to mm/m2/day using
     *     ideal gas equation
     */
    return adt/cmass*8.314*degCtoK(temp)/p*1000.0;
  }
} /* of 'photosynthesis' */
