/***************************************************************************/
/**                                                                       **/
/**               l  e  f  t  m  o  s  t  z  e  r  o  .  c                **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 11.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"
#include "numeric.h"

#define NSEG 20


Real leftmostzero(Real (*fcn)(Real,void *),
                  Real x1,
                  Real x2,
                  void *data,
                  Real xacc,
                  Real yacc,
                  int maxiter)
{
  Real dx,xmid,swap;
  if(x2<x1)
  {
    swap=x1;
    x1=x2;
    x2=swap;
  }

  dx=(x2-x1)/NSEG;
  if((*fcn)(x1,data)<0)
    for(xmid=x1+dx;(*fcn)(xmid,data)<0 && xmid<=x2-dx;xmid+=dx);  
  else
    for(xmid=x1+dx;(*fcn)(xmid,data)>0 && xmid<=x2-dx;xmid+=dx);
  return bisect(fcn,xmid-dx,xmid,data,xacc,yacc,maxiter);
}  /* of 'leftmostzero' */ 
