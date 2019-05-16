/***************************************************************************/
/**                                                                       **/
/**                   b  i  s  e  c  t  .  c                              **/
/**                                                                       **/
/**     Finds a zero of a function using the bisectioning algorithm       **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 11.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include <math.h>
#include "types.h"   /* Definition of datatype Real  */
#include "numeric.h"

Real bisect(Real (*fcn)(Real,void *), /* function */
            Real xlow,  /* lower bound of interval */
            Real xhigh, /* upper bound of interval */
            void *data, /* pointer to additional data for function */
            Real xacc,  /* accuracy in x */
            Real yacc,  /* accuracy in y */
            int maxit   /* maximum number of iterations */
           )            /* returns position of zero of function */
{
  int i;
  Real ylow,yhigh,ymid,xmid;
  ylow=(*fcn)(xlow,data); 
  yhigh=(*fcn)(xhigh,data);
/*#ifdef SAFE
  if(ylow*yhigh>0)
    fail("No zero in [%g,%g]",xlow,xhigh);
#endif*/
  for(i=0;i<maxit;i++)
  {
    xmid=(xlow+xhigh)*0.5;
    if(xhigh-xlow<xacc)
      return xmid;
    ymid=(*fcn)(xmid,data);
    if(fabs(ymid)<yacc)
      return xmid;
    if(ylow*ymid<=0)
    {
      xhigh=xmid;
      yhigh=ymid;
    }
    else
    {
      xlow=xmid;
      ylow=ymid;
    } 
  } /* of for */
  return xmid;
} /* of 'bisect' */
