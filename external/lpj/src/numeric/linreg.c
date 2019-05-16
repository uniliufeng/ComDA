/***************************************************************************/
/**                                                                       **/
/**                      l  i  n  r  e  g  .  c                           **/
/**                                                                       **/
/**     Calculate coefficients of linear regression                       **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 05.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"
#include "numeric.h"

void linreg(Real *a,Real *b,const Real y[],int n)
{
  Real sx,sy,sxx,sxy,ri,delta;
  int i;
  sx=(n*(n-1))/2; /* sum of i from 0 .. n-1 */
  sxx=(n*(n-1)*(2*n-1))/6; /* sum of i^2 from 0 .. n-1 */
  sy=sxy=ri=0;
  for(i=0;i<n;i++)
  {
    ri++;
    sy+=y[i];
    sxy+=y[i]*ri;
  }
  delta=1/(ri*sxx-sx*sx);
  *a=(sxx*sy-sx*sxy)*delta;
  *b=(ri*sxy-sx*sy)*delta;
} /* of 'linreg' */
