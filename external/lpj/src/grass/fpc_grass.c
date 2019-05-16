/***************************************************************************/
/**                                                                       **/
/**             f  p  c  _  g  r  a  s  s  .  c                           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 24.08.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

Real fpc_grass(Pft *pft)
{
  Real fpc_old;
  fpc_old=pft->fpc;
  pft->fpc=pft->nind*(1.0-exp(-0.5*lai_grass(pft)));
  return (pft->fpc<fpc_old) ? 0 : pft->fpc-fpc_old;
} /* 'fpc_grass' */
