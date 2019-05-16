/***************************************************************************/
/**                                                                       **/
/**             p  h  e  n  o  l  o  g  y  _  t  r  e  e  .  c            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"

#define COLDEST_DAY_NHEMISPHERE 14
#define COLDEST_DAY_SHEMISPHERE 195
#define APHEN_MAX 210

Real phenology_tree(Pft *pft,Real temp,Real lat,int day)
{
  Pfttree *tree;
  Real dtemp,gddtw,phen;
  
  tree=pft->data;
  dtemp=temp - getpftpar(pft,gddbase);
  gddtw=temp - getpftpar(pft,twmax);
  if(dtemp>0.0)
    pft->gdd+=dtemp;
  tree->gddtw+= (gddtw>0.0) ? gddtw : 0.0;  
  switch(getpftpar(pft,phenology))
  {
    case SUMMERGREEN: 
      if(pft->aphen<APHEN_MAX)
      {  
        pft->phen=pft->gdd*getpftpar(pft,ramp);
       if(pft->phen>1)
         pft->phen=1;
      }
      else
        pft->phen=0.0;
      break;
    case RAINGREEN: 
      if(pft->wscal<getpftpar(pft,minwscal)) 
        pft->phen=0.0;
      else
      {
        pft->phen=1.0;
        tree->aphen_raingreen++;
      }
      break;
    case EVERGREEN:
      pft->phen=1;
      break;
  } /* of 'switch' */
  if ((lat>=0.0 && day==COLDEST_DAY_NHEMISPHERE) ||
      (lat<0.0 && day==COLDEST_DAY_SHEMISPHERE)) 
    pft->aphen=pft->gdd=0.0;
  pft->aphen+=pft->phen;
  phen=pft->phen;
  return phen;
} /* of 'phenology_tree' */
