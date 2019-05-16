/***************************************************************************/
/**                                                                       **/
/**                     u  p  d  a  t  e  g  d  d  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 26.05.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void updategdd(Real gdd[],
               const Pftpar pftpar[], /* PFT parameter vector of length npft */
               int npft,              /* Number of PFTs */
               Real temp              /* temperature (deg C) */
              )
{
  int p;
  Real dtemp;
  for(p=0;p<npft;p++)
  {
    dtemp=temp-pftpar[p].gddbase;
    if(dtemp>0)
      gdd[p]+=dtemp;
  }
} /* of 'updategdd' */
