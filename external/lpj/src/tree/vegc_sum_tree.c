/***************************************************************************/
/**                                                                       **/
/**                 v  e  g  c  _  s  u  m  _  t  r  e  e  .  c           **/
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

Real vegc_sum_tree(const Pft *pft)
{
  const Pfttree *tree;
  tree=pft->data;
  return (phys_sum_tree(tree->ind)-tree->ind.debt)*pft->nind;
} /* of 'vegc_sum_tree' */
