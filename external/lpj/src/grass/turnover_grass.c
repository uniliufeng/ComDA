/***************************************************************************/
/**                                                                       **/
/**           t  u  r  n  o  v  e  r  _  g  r  a  s  s  .  c              **/
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
#include "grass.h"

/*
 *
 *  Function turnover
 *
 *  Turnover of PFT-specific fraction from each living C pool
 *  Leaf and root C transferred to litter, sapwood C to heartwood
 *
 */

Real turnover_grass(Litter *litter,  /* Litter */
              Pft *pft /* Pointer to pft */
             )         /* returns turnover */
{
  Pftgrass *grass;
  const Pftgrasspar *grasspar;
  Grassphys gturn;
  grass=pft->data;
  grasspar=getpftpar(pft,data);
  gturn.leaf=grass->ind.leaf/grasspar->turnover.leaf; 
  gturn.root=grass->ind.root/grasspar->turnover.root;
  grass->ind.leaf-= gturn.leaf;
  grass->ind.root-= gturn.root;
  litter->ag_grass+=gturn.leaf*pft->nind;
  litter->bg+=gturn.root*pft->nind;
  return gturn.leaf+gturn.root;
} /* of 'turnover_grass' */
