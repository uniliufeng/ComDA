/***************************************************************************/
/**                                                                       **/
/**                    o  u  t  p  u  t  .  c                             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

FILE **fopenoutput(Config config)
{
  int i;
  FILE **output;
  output=newvec(FILE *,config.n_out);
  for(i=0;i<config.n_out;i++)
   if((output[i]=fopen(config.out_filename[i],"wb"))==NULL)
     printfopenerr("fopenoutput",config.out_filename[i]);
  return output;
} /* of 'fopenoutput' */

void fcloseoutput(FILE **output,int n_out)
{
  int i;
  for(i=0;i<n_out;i++)
    if(output[i]!=NULL)
      fclose(output[i]);
  free(output);
} /* of 'fcloseoutput' */

static int writevec(FILE *file,const Real data[],int n,Real frac)
{
  float *vec;
  int i,rc;
  if(file==NULL)
    return 0;
  vec=newvec(float,n);
  for(i=0;i<n;i++)
    vec[i]=(float)(data[i]*frac);
  rc=fwrite(vec,sizeof(float)*n,1,file);
  free(vec);
  return rc;
} /* of 'writevec' */

static int writedata(FILE *file,Real data)
{
  float s;
  if(file==NULL)
   return 0;
  s=(float)(data);
  return fwrite(&s,sizeof(s),1,file);
} /* of 'writedata' */

static int asiwritedata(FILE *file,Real data)
{
  float s;
  if(file==NULL)
   return 0;
  s=(float)(data);
  return fprintf(file,"%f\n",s);
}

int fwriteoutput_annual(FILE **files,const Cell *cell,int n_out,int npft)
{
  
  int count,s,p;
  int i;  // add by wxf for cycle
  Stand *stand;
  Pft *pft;
  float n;
  Real *nfpc;
  Real soilc,litc,vegc,frac;

  vegc=litc=soilc=0.0;
  nfpc=newvec(Real,npft);
  frac=(cell->standlist->n==0) ? 0 : 1.0/(Real)cell->standlist->n;
  n=cell->standlist->n;
  count=0;
 // if(files[FPC]!=NULL)
 //   count+=fwrite(&n,sizeof(n),1,files[FPC]);
 
  if(files[FPC]!=NULL)  //add by wxf
  {	 
 //   count+=fprintf(files[FPC],"%d\n",n);
  }  //add by wxf

  foreachstand(stand,s,cell->standlist)
  {
    for(p=0;p<npft;p++)
      nfpc[p]=0;
    foreachpft(pft,p,stand->pftlist)
    {
       vegc+=vegc_sum(pft);
       nfpc[getpftpar(pft,id)]=pft->fpc;
    }
	for (i=0;i<npft;i++)  // add by wxf
	{
	  	fprintf(files[FPC],"%f,",nfpc[i]);///////////
	}
      	fprintf(files[FPC],"\n"); // add by wxf


    //writevec(files[FPC],nfpc,npft,100.0);
    soilc+=stand->soil.cpool.slow+stand->soil.cpool.fast+stand->soil.litter.bg;
    litc+=stand->soil.litter.ag_tree+stand->soil.litter.ag_grass;
  }
  free(nfpc);
  count+=asiwritedata(files[FIREC],cell->output.firec*frac);
  count+=asiwritedata(files[FIREF],cell->output.firef*frac);
  count+=asiwritedata(files[VEGC],vegc*frac);
  count+=asiwritedata(files[SOILC],soilc*frac);
  count+=asiwritedata(files[LITC],litc*frac);
  count+=asiwritedata(files[FLUX_ESTAB],cell->output.flux_estab*frac);
  return count;
} /* of 'fwriteoutput_annual' */

int fwriteoutput_monthly(FILE **files,const Cell *cell,int n_out,int npft)
{
  
  int count,s,p;
  Stand *stand;
  Pft *pft;
  Real n,frac;
  frac=(cell->standlist->n==0) ? 0 : 1.0/(Real)cell->standlist->n;
  n=cell->standlist->n;
  count=n;
  
  count+=asiwritedata(files[MNPP],cell->output.mnpp*frac);
  count+=asiwritedata(files[MRH],cell->output.mrh*frac);
  count+=asiwritedata(files[MRUNOFF],cell->output.mrunoff*frac);
  count+=asiwritedata(files[MTRANSP],cell->output.mtransp*frac);
  count+=asiwritedata(files[MEVAP],cell->output.mevap*frac);
  count+=asiwritedata(files[MINTERC],cell->output.minterc*frac);
  count+=asiwritedata(files[MSWC1],cell->output.mswc[0]*frac);
  count+=asiwritedata(files[MSWC2],cell->output.mswc[1]*frac);
  return count;
} /* of 'fwriteoutput_monthly' */

void initoutput_annual(Output *output)
{
    int l;
    
    output->mnpp=output->mrh=output->mtransp=output->mrunoff=
    output->mevap=output->minterc=output->flux_estab=0;
    for(l=0;l<NSOILLAYER;l++)
      output->mswc[l]=0;
    output->firec=output->firef=0;
} /* of 'initoutput_annual' */
         
void initoutput_monthly(Output *output)
{
    int l;
    
    output->mnpp=output->mrh=output->mtransp=output->mrunoff=
    output->mevap=output->minterc=0;
    for(l=0;l<NSOILLAYER;l++)
      output->mswc[l]=0;
} /* of 'initoutput_monthly' */
         
