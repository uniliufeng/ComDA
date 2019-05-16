/***************************************************************************/
/**                                                                       **/
/**     l  i  t  t  e  r  _  u  p  d  a  t  e  _  g  r  a  s  s  .  c     **/
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

void litter_update_grass(Litter *litter,const Pft *pft,Real frac)
{
  const Pftgrass *grass;
  grass=pft->data;
  litter->ag_grass+=grass->ind.leaf*frac;
  litter->bg+=grass->ind.root*frac;
#ifdef DEBUG
  printf("ag_grass=%.2f bg=%.2f \n",litter->ag_grass,litter->bg);
#endif
} /* of 'litter_update_grass' */
