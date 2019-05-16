/***************************************************************************/
/**                                                                       **/
/**     l  i  t  t  e  r  _  u  p  d  a  t  e  _  t  r   e  e  .  c       **/
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
#include "tree.h"

void litter_update_tree(Litter *litter,const Pft *pft,Real frac)
{
  const Pfttree *tree;
  tree=pft->data;
  litter->ag_tree+=(tree->ind.leaf+tree->ind.sapwood+tree->ind.heartwood-
                    tree->ind.debt)*frac;
  litter->bg+=tree->ind.root*frac;
} /* of 'litter_update_tree' */
