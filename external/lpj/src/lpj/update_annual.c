/***************************************************************************/
/**                                                                       **/
/**               u  p  d  a  t  e  _  a  n  n  u  a  l  .  c             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.03.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "tree.h"


#define APREC_MIN  100	/* minimum precipitation to establish*/

Real update_annual(FILE **output,
                   Cell *cell,
                   const Pftpar pftpar[],
                   int npft,
                   int ntypes,
                   int year,
                   Config config
                  )
{
  int s,p,m;
  Pft *pft;
  Real *fpc_inc,bm_inc;
  Real totc1,delta_totc;
  Real gdd20,turnover_ind,fire_flux,acflux_estab,fire_frac;
  Stand *stand;
  annual_climbuf(&cell->climbuf);
  
  delta_totc=0.0;
  foreachstand(stand,s,cell->standlist)
  {
    totc1=acflux_estab=0;
    fpc_inc=newvec(Real,getnpft(stand->pftlist));

    foreachpft(pft,p,stand->pftlist)
    {
      bm_inc=pft->bm_inc;
      /*if(isneg(pft,bm_inc))
      { 
        litter_update(&stand->soil.litter,pft,pft->nind);
        delpft(stand->pftlist,p);
        p--;  
        continue;
      }*/

      reproduction(&stand->soil.litter,pft);
      turnover_ind=turnover(&stand->soil.litter,pft);

#ifdef DEBUG
      printf("PFT:%s fpc_inc=%g fpc=%g\n",pft->par->name,fpc_inc[p],pft->fpc);
      printf("PFT:%s %d  bm_inc=%g vegc=%g soil=%g\n",pft->par->name,year,pft->bm_inc,vegc_sum(pft),
      stand->soil.cpool.slow+stand->soil.cpool.fast+stand->soil.litter.bg+stand->soil.litter.ag_tree+stand->soil.litter.ag_grass);
#endif
     
     allocation(&stand->soil.litter,pft,fpc_inc+p);

      bm_inc=0.0;
      if(isneg(pft,bm_inc))
      { 
        fpc_inc[p]=fpc_inc[getnpft(stand->pftlist)-1];
        litter_update(&stand->soil.litter,pft,pft->nind);
        delpft(stand->pftlist,p);
        p--;  
        continue;
      } 
      mortality(&stand->soil.litter,pft,turnover_ind,
                cell->climbuf.temp_max,&fpc_inc[p]);
      gdd20=updatebuffer(pft->gddbuf,cell->gdd[getpftpar(pft,id)]);
      if(isneg(pft,bm_inc))
      { 
        fpc_inc[p]=fpc_inc[getnpft(stand->pftlist)-1];
        litter_update(&stand->soil.litter,pft,pft->nind);
        delpft(stand->pftlist,p);
        p--;  
        continue;
      }
      if(killpft(&stand->soil.litter,pft,&cell->climbuf))
      {
        fpc_inc[p]=fpc_inc[getnpft(stand->pftlist)-1];
        litter_update(&stand->soil.litter,pft,pft->nind);
        delpft(stand->pftlist,p);
        p--;  /* adjust loop index */
#ifdef DEBUG
        printf("killpft\n");
#endif
        continue;
      }
    } /* of foreachpft */
#ifdef DEBUG
    printf("Number of updated pft: %d\n",stand->pftlist->n);
#endif

    light(&stand->soil.litter,stand->pftlist,ntypes,fpc_inc);

    if(iffire)
      fire_flux=firepft(&stand->soil.litter,stand->pftlist,
                         fire_prob(&stand->soil.litter,stand->fire_sum,&fire_frac));

    if (cell->aprec>=APREC_MIN)
     acflux_estab=establishmentpft(&stand->soil.litter,stand->pftlist,pftpar,
                               npft,ntypes,cell->gdd,&cell->climbuf);

    cell->output.flux_estab+=acflux_estab;
    cell->output.firef+=fire_frac;
    cell->output.firec+=fire_flux;
   
    foreachpft(pft,p,stand->pftlist)
     totc1+=vegc_sum(pft);
    totc1+=stand->soil.cpool.slow+stand->soil.cpool.fast+
           stand->soil.litter.bg+stand->soil.litter.ag_tree+stand->soil.litter.ag_grass;
    delta_totc+=(totc1-stand->totc2)*stand->frac;    
    stand->totc2=0;
    foreachpft(pft,p,stand->pftlist)
     stand->totc2+=vegc_sum(pft);
    stand->totc2+=stand->soil.cpool.slow+stand->soil.cpool.fast+
                  stand->soil.litter.bg+stand->soil.litter.ag_tree+stand->soil.litter.ag_grass;

//    if(year>=config.firstyear) 
      fwriteoutput_annual(output,cell,config.n_out,npft);
    free(fpc_inc);
  } /* of foreachstand */
  return delta_totc;
} /* of 'update_annual' */
