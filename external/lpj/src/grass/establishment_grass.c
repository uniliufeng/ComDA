/***************************************************************************/
/**                                                                       **/
/**    e  s  t  a  b  l  i  s  h  m  e  n  t  _  g  r  a  s  s  .  c      **/
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
#include "grass.h"

Real establishment_grass(Pft *pft,Real fpc_total,
                         Real fpc_type,int n_est) 
{

  Pftgrass *grass;
  Pftgrasspar *grasspar;
  Real est_pft,acflux_est;
  if(n_est>0)
  {
    grass=pft->data;
    grasspar=getpftpar(pft,data);
    est_pft=(1.0-fpc_total)/(Real)n_est;
          /* Account for flux from atmosphere to grass regeneration*/

    acflux_est=phys_sum_grass(grasspar->sapl)*est_pft;

          /* Add regeneration biomass to overall biomass*/

    grass->ind.leaf+=grasspar->sapl.leaf*est_pft;
    grass->ind.root+=grasspar->sapl.root*est_pft;
  }
  else
    acflux_est=0;
  fpc_grass(pft);
  return acflux_est;
} /* of 'establishment_grass' */
