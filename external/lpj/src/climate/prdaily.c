/***************************************************************************/
/**                                                                       **/
/**                   p  r  d  a  i  l  y  .  c                           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 04.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

/*
 *  PRDAILY
 *  Distribution of monthly precipitation totals to quasi-daily values
 *  
 */

#define c1 1.0	/*normalising coefficient for exponential distribution*/
#define c2 1.2	/*power for exponential distribution*/

void prdaily(Climbuf *climbuf,  /* Climate buffer */
             int ndaymonth,     /* number of days in month */
             Real mval,	        /* total rainfall (mm) for month */
             Real mval_wet 	/* expected number of rain days for month */
            )
{


  int d;
  Real prob_rain; 	/* daily probability of rain for this month */
  Real mprec; 		/* average rainfall per rain day for this month */
  Real mprec_sum; 	/* cumulative sum of rainfall for this month */
  Real prob;


  if (mval<1.0) 
    for (d=1;d<=ndaymonth;d++) 
      climbuf->dval_prec[d]=0.0;
  else 
  {
    mprec_sum=0.0;
    if (mval_wet<1.0)
      mval_wet=1.0;
    prob_rain=mval_wet/ndaymonth;
    mprec=mval/mval_wet;

    while (mprec_sum<1.0)
    {
      for (d=1;d<=ndaymonth;d++) 
      {
      
/*----------Transitional probabilities (Geng et al 1986)--------*/

        prob=(climbuf->dval_prec[d-1]<0.1) ? 0.75*prob_rain : 
                                               0.25+(0.75*prob_rain);

/*---------Determine wet days randomly and use Krysanova/Cramer estimates of
	    parameter values (c1,c2) for an exponential distribution---------*/
#ifdef USE_RAND48
        if (drand48()>prob) 
#else
        if (randfrac()>prob) 
#endif
        climbuf->dval_prec[d]=0.0;
        else 
        {
#ifdef USE_RAND48
 	  climbuf->dval_prec[d]=pow(-log(drand48()),c2)*mprec*c1;
#else
 	  climbuf->dval_prec[d]=pow(-log(randfrac()),c2)*mprec*c1;
#endif
          if (climbuf->dval_prec[d]<0.1) 
            climbuf->dval_prec[d]=0.0;
          mprec_sum+=climbuf->dval_prec[d];
        }
      } /* of 'for(d=1;...)' */
     } 

/*-------- Normalise generated precipitation by prescribed monthly totals----*/

      if (mprec_sum>1.0) 
      {
        for (d=1;d<=ndaymonth;d++) 
        {
          climbuf->dval_prec[d]*=mval/mprec_sum;
          if (climbuf->dval_prec[d]<0.1) 
            climbuf->dval_prec[d]=0.0;
        }
      }
  }
  climbuf->dval_prec[0]=climbuf->dval_prec[ndaymonth];
} /* of 'prdaily' */
