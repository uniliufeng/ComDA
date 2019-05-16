/***************************************************************************/
/**                                                                       **/
/**                  t  u  r  n  o  v  e  r  _  t  r  e  e  .  c          **/
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

/*
 *
 *  Function turnover
 *
 *  Turnover of PFT-specific fraction from each living C pool
 *  Leaf and root C transferred to litter, sapwood C to heartwood
 *
 */

Real turnover_tree(Litter *litter,  /* Litter           */
                   Pft *pft         /* Pointer to pft   */
                  )                 /* returns turnover */
{
  Pfttree *tree;
  const Pfttreepar *treepar;
  Treephys turn;
  tree=pft->data;
  treepar=getpftpar(pft,data);
  turn.leaf=tree->ind.leaf/treepar->turnover.leaf; 
  turn.sapwood=tree->ind.sapwood/treepar->turnover.sapwood; 
  turn.root=tree->ind.root/treepar->turnover.root;
  tree->ind.leaf-= turn.leaf;
  tree->ind.sapwood-= turn.sapwood;
  tree->ind.root-= turn.root;
  tree->ind.heartwood+=turn.sapwood;
  litter->ag_tree+=turn.leaf*pft->nind;
  litter->bg+=turn.root*pft->nind;
#ifdef DEBUG2
  printf("turn: %g %g %g %g %g %g %g %g %g\n",tree->ind.leaf,tree->ind.root,
         tree->ind.sapwood,tree->ind.heartwood,tree->ind.debt,pft->fpc,pft->nind,tree->height,
         tree->crownarea,pft->bm_inc);
#endif
  return turn.leaf+turn.sapwood+turn.root;
} /* of 'turnover_tree' */
