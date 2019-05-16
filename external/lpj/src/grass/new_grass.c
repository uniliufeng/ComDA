/***************************************************************************/
/**                                                                       **/
/**                       n  e  w  _  g  r  a  s  s  .  c                 **/
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

void new_grass(Pft *pft /* Parameter of pft */
              )        
{
  Pftgrass *grass;
  grass=new(Pftgrass);
  pft->data=grass;
  pft->nind=1;
  grass->ind.leaf=grass->ind.root=0;
} /* of 'new_grass' */
