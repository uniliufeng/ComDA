/***************************************************************************/
/**                                                                       **/
/**                 f  w  r  i  t  e  _  t  r  e  e  .  c                 **/
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

Bool fwrite_tree(FILE *file,const Pft *pft)
{
  const Pfttree *tree;
  tree=pft->data;
  return fwrite(tree,sizeof(Pfttree),1,file)!=1;
} /* of 'fwrite_tree' */
