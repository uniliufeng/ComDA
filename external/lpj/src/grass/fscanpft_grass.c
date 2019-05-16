/***************************************************************************/
/**                                                                       **/
/**               f  s  c  a  n  p  f  t  _  g  r  a  s  s  .  c          **/
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
#include "grass.h"

#define fscangrassphys2(file,var,fcn,name) \
  if(fscangrassphys(file,var))\
  { \
    fprintf(stderr,"Error reading limit '%s' in '%s'.\n",name,fcn); \
    return TRUE; \
  }


static Bool fscangrassphys(FILE *file,Grassphys *phys)
{
  double leaf,root;
  if(fscanf(file,"%lg %lg",&leaf,&root)!=2)
    return TRUE;
  if(leaf<=0 ||  root<=0)
    return TRUE;
  phys->leaf=leaf;
  phys->root=root;
  return FALSE;
} /* of 'fscangrassphys' */

Bool fscanpft_grass(FILE *file,          /* file pointer */
                    Pftpar *pft,         /* Pointer to Pftpar array */
                    const char *filename /* filename */
                   )                     /* returns FALSE for success  */
{
  Pftgrasspar *grass;
  pft->newpft=new_grass;
  pft->turnover=turnover_grass;
  pft->npp=npp_grass;
  pft->fpc=fpc_grass;
  pft->leaf_phenology=phenology_grass;
  pft->fwrite=fwrite_grass;
  pft->fprint=fprint_grass;
  pft->fread=fread_grass;
  pft->litter_update=litter_update_grass;
  pft->allocation=allocation_grass;
  pft->isneg=isneg_grass;
  pft->establishment=establishment_grass;
  pft->reproduction=reproduction_grass;
  pft->lai=lai_grass;
  pft->init=init_grass;
  pft->free=free_grass;
  pft->vegc_sum=vegc_sum_grass;
  pft->light=light_grass;
  grass=new(Pftgrasspar);
  pft->data=grass;
  fscangrassphys2(file,&grass->turnover,filename,"turnover");
  fscangrassphys2(file,&grass->cn_ratio,filename,"cn_ratio");
  grass->sapl.leaf=pft->lai_sapl/pft->sla;
  grass->sapl.root=(1.0/pft->lmro_ratio)*grass->sapl.leaf;
  return FALSE;
} /* of 'fscanpft_grass' */
