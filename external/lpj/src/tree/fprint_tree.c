/***************************************************************************/
/**                                                                       **/
/**                 f  p  r  i  n  t  _  t  r  e  e  .  c                 **/
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


void fprint_tree(FILE *file,const Pft *pft)
{
  Pfttree *tree;
  tree=pft->data;
  fprintf(file,"Cmass: ");
  fprinttreephys2(file,tree->ind,pft->nind);
  fprintf(file,"\n");
} /* of 'fprint_tree' */
