/***************************************************************************/
/**                                                                       **/
/**               f  r  e  a  d  _  g  r  a  s  s  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"

void fread_grass(FILE *file,Pft *pft,Bool swap)
{
  Pftgrass *grass;
  grass=new(Pftgrass);
  pft->data=grass;
  freadreal((Real *)&grass->ind,sizeof(Grassphys)/sizeof(Real),swap,file);
  fpc_grass(pft);
} /* of 'fread_grass' */
