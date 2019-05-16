/***************************************************************************/
/**                                                                       **/
/**                 f  p  r  i  n  t  _  g  r  a  s  s  .  c              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"


void fprint_grass(FILE *file,const Pft *pft)
{
  Pftgrass *grass;
  grass=pft->data;
  fprintf(file,"Cmass: ");
  fprintgrassphys(file,grass->ind);
  fprintf(file,"\n");
} /* of 'fprint_grass' */
