/***************************************************************************/
/**                                                                       **/
/**               f  w  r  i  t  e  s  o  i  l  .  c                      **/
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

Bool fwritesoil(FILE *file,const Soil *soil,Bool full)
{
  fwrite1(&soil->cpool,sizeof(Pool),file);
  fwrite1(&soil->litter,sizeof(Litter),file);
  fwriten(soil->w,sizeof(Real),NSOILLAYER,file);
  fwrite1(&soil->w_evap,sizeof(Real),file);
  fwrite1(&soil->snowpack,sizeof(Real),file);
  if(full)
  {
    fwrite1(&soil->alag,sizeof(Real),file);
    fwrite1(&soil->amp,sizeof(Real),file);
    fwrite1(&soil->meanw1,sizeof(Real),file);
  }
  return FALSE;
} /* of 'fwritesoil' */
