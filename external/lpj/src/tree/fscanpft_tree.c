/***************************************************************************/
/**                                                                       **/
/**               f  s  c  a  n  p  f  t  _  t  r  e  e  .  c             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"


#define fscanreal2(file,var,fcn,pft,name) \
  if(fscanreal(file,var)) \
  { \
    fprintf(stderr,"Error reading '%s' of PFT '%s' in '%s'.\n",name,pft,fcn); \
    return TRUE; \
  }
#define fscanint(file,var,fcn,pft,name) \
  if(fscanf(file,"%d",var)!=1) \
  { \
    fprintf(stderr,"Error reading '%s' of PFT '%s' in '%s'.\n",name,pft,fcn); \
    return TRUE; \
  }
#define fscantreephys2(file,var,fcn,name)\
  if(fscantreephys(file,var))\
  {\
    fprintf(stderr,"Error reading treephys '%s' in '%s'.\n",name,fcn); \
    return TRUE; \
  }

static Bool fscantreephys(FILE *file,Treephys *phys)
{
  double leaf,sapwood,root;
  if(fscanf(file,"%lg %lg %lg",&leaf,&sapwood,&root)!=3)
    return TRUE;
  if(leaf<=0 || sapwood<=0 || root<=0)
    return TRUE;
  phys->leaf=leaf;
  phys->sapwood=sapwood;
  phys->root=root;
  return FALSE;
} /* of 'fscantreephys' */

Bool fscanpft_tree(FILE *file,           /* file pointer */
                   Pftpar *pft,          /* Pointer to Pftpar array */
                   const char *filename  /* filename */
                  )                      /* returns  FALSE for success */
{
  Real stemdiam,height_sapl,wood_sapl;
  Pfttreepar *tree;
  pft->newpft=new_tree;
  pft->turnover=turnover_tree;
  pft->npp=npp_tree;
  pft->fpc=fpc_tree;
  pft->leaf_phenology=phenology_tree;
  pft->fwrite=fwrite_tree;
  pft->fprint=fprint_tree;
  pft->fread=fread_tree;
  pft->litter_update=litter_update_tree;
  pft->allocation=allocation_tree;
  pft->isneg=isneg_tree;
  pft->establishment=establishment_tree;
  pft->reproduction=reproduction_tree;
  pft->mortality=mortality_tree;
  pft->init=init_tree;
  pft->fire=fire_tree;
  pft->lai=lai_tree;
  pft->free=free_tree;
  pft->vegc_sum=vegc_sum_tree;
  pft->light=light_tree;
  pft->adjust=adjust_tree;
  tree=new(Pfttreepar);
  pft->data=tree;
  fscanint(file,&tree->leaftype,filename,pft->name,"leaftype");
  fscantreephys2(file,&tree->turnover,filename,"turnover");
  fscantreephys2(file,&tree->cn_ratio,filename,"cn_ratio");
  fscanreal2(file,&tree->crownarea_max,filename,pft->name,"crownarea_max");
  fscanreal2(file,&wood_sapl,filename,pft->name,"wood_sapl");
  tree->sapl.leaf=pow(pft->lai_sapl*allom1*pow(wood_sapl,reinickerp)*
                  pow(4*pft->sla/M_PI/k_latosa,reinickerp*0.5)/pft->sla,
                  2/(2-reinickerp));
  stemdiam=wood_sapl*sqrt(4*tree->sapl.leaf*pft->sla/M_PI/k_latosa);
  height_sapl=allom2*pow(stemdiam,allom3);
  tree->sapl.sapwood=wooddens*height_sapl*tree->sapl.leaf*pft->sla/k_latosa;
  tree->sapl.heartwood=(wood_sapl-1)*tree->sapl.sapwood;
  tree->sapl.root=(1.0/pft->lmro_ratio)*tree->sapl.leaf;
  return FALSE;
} /* of 'fscanpft_tree' */
