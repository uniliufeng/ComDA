/***************************************************************************/
/**                                                                       **/
/**              l  i  g  h  t  _  t  r  e  e  .  c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

void light_tree(Litter *litter,Pft *pft,Real excess)
{
  Pfttree *tree;
  tree=pft->data;

  Real nind_kill;  /* reduction in individual density to reduce tree FPC to
                      permitted maximum (indiv/m2)*/
#ifdef DEBUG2*/
  printf("light: %g %g %s %g\n",pft->fpc,pft->nind,pft->par->name,excess);
#endif*/
  nind_kill=(excess<1e-20) ? 0 : pft->nind*(excess/pft->fpc); 
  pft->nind-=nind_kill;
  litter_update_tree(litter,pft,nind_kill);
  fpc_tree(pft);

} /* of 'light_tree' */
