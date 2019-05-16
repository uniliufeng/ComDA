/***************************************************************************/
/**                                                                       **/
/**               i  s  n  e  g  _  t  r  e  e  .  c                      **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.04.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

Bool isneg_tree(const Pft *pft,Real bm_inc)
{
  Pfttree *tree;
  tree=pft->data;
#ifdef DEBUG
  printf("isneg: %g %g %g %g %g %g %g %g\n",tree->ind.leaf*pft->nind,tree->ind.root*pft->nind,
         tree->ind.sapwood*pft->nind,tree->ind.heartwood*pft->nind,tree->ind.debt*pft->nind,pft->fpc,pft->nind,pft->bm_inc);
#endif
  return ((tree->ind.leaf+tree->ind.root+tree->ind.sapwood+tree->ind.heartwood-tree->ind.debt)<0.0 
         || tree->ind.root<0.0 || tree->ind.leaf<0.0 || tree->ind.sapwood<0.0 || tree->ind.heartwood<0.0
         ||pft->fpc<=1e-20 || pft->nind<=1e-20);
} /* of 'isneg_tree' */
