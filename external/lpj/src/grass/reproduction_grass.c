/***************************************************************************/
/**                                                                       **/
/**          r  e  p  r  o  d  u  c  t  i  o  n  _  g  r  a  s  s  . c    **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 28.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

#define reprod_cost 0.1

void reproduction_grass(Litter *litter,Pft *pft)
{
  Real reprod;
  if(pft->bm_inc>0)
  {
    reprod=pft->bm_inc*reprod_cost;
    litter->ag_grass+=reprod;
    pft->bm_inc-=reprod;
  }
} /* of 'reproduction_grass' */
