/***************************************************************************/
/**                                                                       **/
/**                 f  w  r  i  t  e  _  g  r  a  s  s  .  c              **/
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

Bool fwrite_grass(FILE *file,const Pft *pft)
{
  const Pftgrass *grass;
  grass=pft->data;
  return fwrite(grass,sizeof(Pftgrass),1,file)!=1;
} /* of 'fwrite_grass' */
