/***************************************************************************/
/**                                                                       **/
/**               e  s  t  a  b  l  i  s  h  m  e  n  t  .  c             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"
#include "grass.h"

Real establishment(Litter *litter,Pftlist pftlist,const Pftpar *pftpar,
                   int npft,int ntypes,const Real gdd[],const Climbuf *climbuf) 
{

  /*
   *  DESCRIPTION
   *
   *  Establishment in population (standard LPJ) mode.
   *  Simulates population increase through establishment each simulation 
   *  year for trees and grasses and introduces new PFTs if climate conditions
   *  become suitable.
   *  This function assumes each Individual object represents the average
   *  individual of a PFT population, and that there is (at most) one 
   *  individual object per PFT per modelled area (stand).
   *
   */

  Real est_pft;	/* establishment rate for a particular PFT on modelled area 
                   basis (for trees, indiv/m2; for grasses, fraction of 
                   modelled area colonised establishment rate for a particular
                   PFT on modelled area basis (for trees, indiv/m2; for 
                   grasses, fraction of modelled area colonised) */
  Real acflux_est;
  Real fpc_total,*fpc_type;
  int *n_est; 
  Bool *present;
  int p,i;
  Pft *pft;
  present=newvec(Bool,npft);
  fpc_type=newvec(Real,ntypes);
  n_est=newvec(int,ntypes);
  acflux_est=0.0;
  for(p=0;p<npft;p++)
    present[p]=FALSE;
  foreachpft(pft,p,pftlist)
    present[pft->par->id]=TRUE;
#ifdef DEBUG
  printf("Number of pfts: %d\n",pftlist->n);
  for(p=0;p<npft;p++)
    printf("%s ",present[p] ? "true" : "false");
  printf("\n");
#endif
  fpc_total=fpc_sum(fpc_type,ntypes,pftlist);
  for(i=0;i<ntypes;i++)
    n_est[i]=0;
  for(p=0;p<npft;p++)
    if(establish(gdd[p],pftpar+p,climbuf))
    {
      if(!present[p])
        addpft(pftlist,pftpar+p);
      n_est[pftpar[p].type]++;
    } 
#ifdef DEBUG
  printf("new npft=%d\n",pftlist->n);
#endif
  foreachpft(pft,p,pftlist)
    if(establish(gdd[pft->par->id],pftpar+pft->par->id,climbuf))
      acflux_est+=pft->par->establishmentpft(pft,fpc_total,fpc_type[pft->par->type],n_est[pft->par->type]);
  fpc_total=fpc_sum(fpc_type,ntypes,pftlist); 
  foreachpft(pft,p,pftlist)
    if(fpc_type[TREE]>FPC_TREE_MAX) pft->par->adjust(litter,pft,fpc_type[pft->par->type]);
#ifdef DEBUG
  printf("new npft=%d\n",pftlist->n);
#endif
  free(present);
  free(fpc_type);
  free(n_est);
  return acflux_est;
} /* of 'establishment' */
