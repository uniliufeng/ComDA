/***************************************************************************/
/**                                                                       **/
/**                             l  i  g  h  t  .  c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 22.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

void light(Litter *litter,Pftlist pftlist,int ntypes,const Real fpc_inc[])

{
  int p;
  Real ntree;       	/*no of tree PFTs currently present*/
  Real *fpc_total; 	/*total grid FPC for PFTs     */ 
  Real fpc_inc_tree;   	/*this years total FPC increment for tree PFTs*/
  Real excess;         	/*tree FPC or grass cover to be reduced*/
  Bool recalc;
  Real nind_kill;
  Pft *pft;
  Real totc1,totc2, soil1,soil2;
/*    Calculate total woody FPC, FPC increment and grass cover (= crown area)*/
 
  totc1=totc2=soil1=soil2=0.0;
  fpc_inc_tree=0.0;
  ntree=0;
  fpc_total=newvec(Real,ntypes);
  fpc_sum(fpc_total,ntypes,pftlist);
/*  foreachpft(pft,p,pftlist) 
    totc1+=vegc_sum(pft);
  soil1=litter->bg+litter->ag_tree+litter->ag_grass;*/
  foreachpft(pft,p,pftlist) 
    if(istree(pft))
    {
      fpc_inc_tree+=fpc_inc[p];
      ntree++;
    }

  foreachpft(pft,p,pftlist) 
  {
    switch(getpftpar(pft,type))
    {
      case TREE:
        if (fpc_total[TREE]>FPC_TREE_MAX) 
        { 
          excess =(fpc_inc_tree>0.0) ? 
                  (fpc_total[TREE]-FPC_TREE_MAX)*(fpc_inc[p]/fpc_inc_tree) :
                  (fpc_total[TREE]-FPC_TREE_MAX)/ntree;

          /*  Reduce individual density (and thereby gridcell-level biomass)*/
          /*  so that total tree FPC reduced to 'fpc_tree_max'*/
          /* changed by Werner von Bloh to avoid FPE */
          
          pft->par->light(litter,pft,excess);
        }
        break;
                              
      case GRASS:
        if(fpc_total[GRASS]>(1.0-min(fpc_total[TREE],FPC_TREE_MAX))) 
        {

          excess=(min(fpc_total[TREE],FPC_TREE_MAX)+fpc_total[GRASS]-1.0)*
                 (pft->fpc/fpc_total[GRASS]);
          pft->par->light(litter,pft,excess);
        }
        break;
    } /* of 'switch' */
  }  /* of 'foreachpft' */
/*  foreachpft(pft,p,pftlist) 
    totc2+=vegc_sum(pft);
  soil2=litter->bg+litter->ag_tree+litter->ag_grass;
  printf("totc1=%g totc2=%g soil1=%g soil2=%g \n",totc1,totc2,soil1,soil2);*/
  free(fpc_total);
} /* of 'light' */
