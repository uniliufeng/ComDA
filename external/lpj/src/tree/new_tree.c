/***************************************************************************/
/**                                                                       **/
/**                       n  e  w  _  t  r  e  e  .  c                    **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 18.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

void new_tree(Pft *pft /* Parameter of pft */
             )         
{
  Pfttree *tree;
  tree=new(Pfttree);
  pft->data=tree;
/*  tree->leafon=FALSE;*/
  tree->ind.root=tree->ind.sapwood=tree->ind.heartwood=tree->ind.leaf=0.0;
  tree->ind.debt=tree->gddtw=tree->aphen_raingreen=0.0;
} /* of 'new_tree' */
