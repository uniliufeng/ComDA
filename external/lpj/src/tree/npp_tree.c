/***************************************************************************/
/**                                                                       **/
/**                    n  p  p  _  t  r  e  e  .  c                       **/
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

#define k 0.0548 

Real npp_tree(Pft *pft,Real phen,Real gtemp_air,Real gtemp_soil,Real assim)
{
  Pfttree *tree;
  const Pfttreepar *par;
  Treephys resp;
  Real mresp,gresp,npp; 
  tree=pft->data;
  par=pft->par->data;
  resp.sapwood=pft->par->respcoeff*k*tree->ind.sapwood*pft->nind/
               par->cn_ratio.sapwood*gtemp_air;
  resp.root=pft->par->respcoeff*k*tree->ind.root*pft->nind/par->cn_ratio.root*
            gtemp_soil*phen;
  gresp=(assim-resp.sapwood-resp.root)*0.25;
  if (gresp<0.0) gresp=0.0;
  mresp=resp.sapwood+resp.root;
#ifdef DEBUG3
  printf("mresp=%g, gresp=%g\n",mresp,gresp);
#endif
  npp=assim-mresp-gresp;
  pft->bm_inc+=npp;
  return npp;
} /* of 'npp_tree' */
