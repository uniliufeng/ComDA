/***************************************************************************/
/**                                                                       **/
/**      e  s  t  a  b  l  i  s  h  m  e  n  t  _  t  r  e  e  .  c       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 28.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

#define K_EST 0.12		/*maximum overall sapling establishment rate (indiv/m2)*/

Real establishment_tree(Pft *pft,Real fpc_total,
                        Real fpc_type,int n_est) 
{

  Real frac;
  Real nind_old,acflux_est,est_pft;
  Pfttree *tree;
  Pfttreepar *treepar;
  tree=pft->data;
  treepar=getpftpar(pft,data);
  
  if (fpc_type>=FPC_TREE_MAX && n_est<=0.0)
  { 
   allometry_tree(pft);
   return 0.0;
  }
  est_pft=K_EST*(1.0-exp(-5.0*(1.0-fpc_type)))*(1.0-fpc_type)/(Real)n_est;
#ifdef SAFE
  if (est_pft<0.0){
     printf("fpc= %g\n",fpc_type);
     fail("establishment_area: negative establishment rate=%g",est_pft);}
#endif
  if (pft->nind<0.0) pft->nind=0.0;
  nind_old=pft->nind;
  pft->nind+=est_pft;

 /*  Account for flux from the atmosphere to new saplings */

  acflux_est=phys_sum_tree(treepar->sapl)*est_pft;

 /* Adjust average individual C biomass based on average biomass and density of the new saplings*/
  tree->ind.leaf=(tree->ind.leaf*nind_old+treepar->sapl.leaf*est_pft)/pft->nind;
  tree->ind.root=(tree->ind.root*nind_old+treepar->sapl.root*est_pft)/pft->nind;
  tree->ind.sapwood=(tree->ind.sapwood*nind_old+treepar->sapl.sapwood*est_pft)/pft->nind;
  tree->ind.heartwood=(tree->ind.heartwood*nind_old+treepar->sapl.heartwood*est_pft)/pft->nind;
  tree->ind.debt=tree->ind.debt*nind_old/pft->nind;
  allometry_tree(pft);
  fpc_tree(pft);
  return acflux_est;
} /* of 'establishment_tree' */
