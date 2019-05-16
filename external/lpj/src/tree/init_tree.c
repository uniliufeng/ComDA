/***************************************************************************/
/**                                                                       **/
/**             i  n  i  t  _  t  r  e  e  .  c                           **/
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

void init_tree(Pft *pft)
{
  Pfttree *tree;
  tree=pft->data;
  tree->gddtw=tree->aphen_raingreen=0;
/*  tree->leafondays=tree->leafoffdays=0;*/
} /* of 'init_tree' */
