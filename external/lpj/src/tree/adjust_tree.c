/***************************************************************************/
/**                                                                       **/
/**               a  d  j  u  s  t  _  t  r  e  e  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 22.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

void adjust_tree(Litter *litter,Pft *pft,Real fpc_tree) 

{
  Real frac;

  frac=FPC_TREE_MAX/fpc_tree;
  pft->nind*=frac;
  pft->fpc*=frac;
  litter_update_tree(litter,pft,(pft->nind/frac-pft->nind)); 
  
} /* of 'adjust_tree' */
