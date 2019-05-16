/***************************************************************************/
/**                                                                       **/
/**             i  n  i  t  _  a  n  n  u  a  l  .  c                     **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 10.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void init_annual(Cell *cell,int npft)
{
  int s,p;
  Pft *pft;
  Stand *stand;
  initgdd(cell->gdd,npft);
  init_climbuf(&cell->climbuf);
  cell->aprec=0.0;
  foreachstand(stand,s,cell->standlist)
  {
#ifdef DEBUG
    printf("init npft=%d\n",stand->pftlist->n); 
#endif
    stand->fire_sum=0;
    foreachpft(pft,p,stand->pftlist)
      initpft(pft);
  } /* of foreachstand */
} /* of 'init_annual' */
