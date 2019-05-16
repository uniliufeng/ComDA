/***************************************************************************/
/**                                                                       **/
/**             l  a  i  _  g  r  a  s  s  .  c                           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

Real lai_grass(const Pft *pft)
{
  return ((Pftgrass *)pft->data)->ind.leaf*getpftpar(pft,sla);
} /* 'lai_grass' */
