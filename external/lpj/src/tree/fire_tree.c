/***************************************************************************/
/**                                                                       **/
/**                    f  i  r  e  _  t  r  e  e  .  c                    **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 28.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

Real fire_tree(Pft *pft,Real fireprob)
{
  Real disturb,flux;
  Pfttree *tree;
  tree=pft->data;
  disturb=(1-pft->par->resist)*fireprob;
  flux=disturb*pft->nind*(tree->ind.leaf+tree->ind.sapwood+
                          tree->ind.heartwood-tree->ind.debt+tree->ind.root);
  pft->nind*=(1-disturb);
  return flux;
} /* of 'fire_tree' */
