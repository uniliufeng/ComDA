/***************************************************************************/
/**                                                                       **/
/**               f  r  e  a  d  _  t  r  e  e  .  c                      **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

void fread_tree(FILE *file,Pft *pft,Bool swap)
{
  Pfttree *tree;
  tree=new(Pfttree);
  pft->data=tree;
/*  freadint1(&tree->leafondays,swap,file);
  freadint1(&tree->leafoffdays,swap,file);
  freadint1(&tree->leafon,swap,file);*/
  freadreal1(&tree->height,swap,file);
  freadreal1(&tree->crownarea,swap,file);
  freadreal1(&tree->gddtw,swap,file);
  freadreal1(&tree->aphen_raingreen,swap,file);
  freadreal((Real *)&tree->ind,sizeof(Treephys2)/sizeof(Real),swap,file);
  fpc_tree(pft);
} /* of 'fread_tree' */
