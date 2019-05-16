/***************************************************************************/
/**                                                                       **/
/**               f  p  r  i  n  t  s  o  i  l  .  c                      **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 28.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void fprintsoil(FILE *file,const Soil *soil)
{
  int l;
  fprintf(file,"Soil type: %s\n",soil->par->name);
  fprintf(file,"Cpool: ");
  fprintpool(file,soil->cpool);
  fprintf(file,"\nLitter: ");
  fprintlitter(file,soil->litter);
  fprintf(file,"\nWater: ");
  foreachsoillayer(l)
    fprintf(file,"%5.2g ",soil->w[l]);
  fprintf(file,"\nSnowpack: %g\n",soil->snowpack);
} /* of 'fprintsoil' */
