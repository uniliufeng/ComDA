/***************************************************************************/
/**                                                                       **/
/**                      f  p  c   _  t  r  e  e  .  c                    **/
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

Real fpc_tree(Pft *pft)
{
  Pfttree *tree;
  Real fpc_old;
  fpc_old=pft->fpc;
  tree=pft->data;
  pft->fpc=(tree->crownarea>0.0) ? tree->crownarea*pft->nind*
                                   (1.0-exp(-0.5*lai(pft))) : 0;
  return (pft->fpc<fpc_old) ? 0 : pft->fpc-fpc_old;
} /* of 'fpc_tree' */
