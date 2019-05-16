/***************************************************************************/
/**                                                                       **/
/**              f  r  e  a  d  s  o  i  l  .  c                          **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool freadsoil(FILE *file,Soil *soil,const Soilpar *soilpar,Bool swap)
{
  soil->par=soilpar;
  freadreal((Real *)&soil->cpool,sizeof(Pool)/sizeof(Real),swap,file);
  freadreal((Real *)&soil->litter,sizeof(Litter)/sizeof(Real),swap,file);
  freadreal(soil->w,NSOILLAYER,swap,file);
  freadreal1(&soil->w_evap,swap,file);
  freadreal1(&soil->snowpack,swap,file);
  freadreal1(&soil->alag,swap,file);
  freadreal1(&soil->amp,swap,file);
  return (freadreal1(&soil->meanw1,swap,file)!=1);
} /* of 'freadsoil' */
