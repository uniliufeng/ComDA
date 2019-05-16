/***************************************************************************/
/**                                                                       **/
/**               f  s  c  a  n  p  f  t  p  a  r  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 05.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define fscanreal2(file,var,fcn,pft,name) \
  if(fscanreal(file,var)) \
  { \
    fprintf(stderr,"Error reading '%s' of PFT '%s' in '%s'.\n",name,pft,fcn); \
    break; \
  }
#define fscanint(file,var,fcn,pft,name) \
  if(fscanf(file,"%d",var)!=1) \
  { \
    fprintf(stderr,"Error reading '%s' of PFT '%s' in '%s'.\n",name,pft,fcn); \
    break; \
  }
#define fscanlimit2(file,var,fcn,pft,name) \
  if(fscanlimit(file,var)) \
  { \
    fprintf(stderr,"Error reading '%s' of PFT '%s' in '%s'.\n",name,pft,fcn); \
    break; \
  }

static Bool fscanlimit(FILE *file,Limit *limit)
{
  double low,high;
  if(fscanf(file,"%lg %lg",&low,&high)!=2)
    return TRUE;
  limit->low=low;
  limit->high=high;
  return FALSE;
} /* of 'fscanlimit' */

int fscanpftpar(Pftpar **pftpar,     /* Pointer to Pftpar array */
                const char *filename, /* filename */
                const Fscanpftparfcn scanfcn[],
                int ntypes           /* number of PFT types */
               )                     /* returns number of elements in array */
{
  int npft,n,id,l;
  String s;
  char *cmd;
  FILE *file;
  Pftpar *pft;
#ifdef USE_CPP
  cmd=malloc(strlen(filename)+strlen(cpp_cmd)+1);
  strcat(strcpy(cmd,cpp_cmd),filename);
  file=pt_popen(cmd,"r");
  free(cmd);
#else
  file=fopen(filename,"r");
#endif
  if(file==NULL)
  {
    printfopenerr("fscanpftpar",filename);
    return 0;
  }

  if(fscanf(file,"%d",&npft)!=1)
  {
    readinterr(filename,"npft");
    pt_pclose(file);
    return 0;
  }
  *pftpar=newvec(Pftpar,npft);
  for(n=0;n<npft;n++)
  {
    if(fscanf(file,"%d",&id)!=1)
    {
      readinterr(filename,"id");
      break;
    } 
    if(id<0 || id>=npft)
    {
      fprintf(stderr,"Error in '%s': invalid range of 'id'=%d.\n",filename,id);
      break;
    } 
    pft=(*pftpar)+id;
    pft->id=id;
    if(fscanstring(file,s))
    {
      readstringerr(filename,"name");
      break;
    }
    pft->name=strdup(s);
    fscanint(file,&pft->type,filename,pft->name,"type");
    pft->rootdist[BOTTOMLAYER]=1;
    for(l=0;l<BOTTOMLAYER;l++)
    {
      fscanreal2(file,pft->rootdist+l,filename,pft->name,"rootdist");
      pft->rootdist[BOTTOMLAYER]-=pft->rootdist[l];
    }
    fscanreal2(file,&pft->minwscal,filename,pft->name,"minwscal");
    fscanreal2(file,&pft->gmin,filename,pft->name,"gmin");
    fscanreal2(file,&pft->respcoeff,filename,pft->name,"respcoeff");
    fscanreal2(file,&pft->nmax,filename,pft->name,"nmax");
    fscanreal2(file,&pft->resist,filename,pft->name,"resist");
    fscanreal2(file,&pft->longivity,filename,pft->name,"longivity");
    pft->sla=2e-4*exp(6.15-0.46*log(pft->longivity*12));
    fscanreal2(file,&pft->lmro_ratio,filename,pft->name,"lmro_ratio");
    fscanreal2(file,&pft->ramp,filename,pft->name,"ramp");
    pft->ramp=1/pft->ramp;
    fscanreal2(file,&pft->lai_sapl,filename,pft->name,"lai_sapl");
    fscanreal2(file,&pft->gdd5min,filename,pft->name,"gdd5min");
    fscanreal2(file,&pft->twmax,filename,pft->name,"twmax");
    fscanreal2(file,&pft->gddbase,filename,pft->name,"gddbase");
    fscanreal2(file,&pft->min_temprange,filename,pft->name,"min_temprange");
    fscanreal2(file,&pft->emax,filename,pft->name,"emax");
    fscanreal2(file,&pft->intc,filename,pft->name,"intc");
    fscanint(file,&pft->phenology,filename,pft->name,"phenology");
    fscanint(file,&pft->path,filename,pft->name,"path");
    fscanlimit2(file,&pft->temp_co2,filename,pft->name,"temp_co2");
    fscanlimit2(file,&pft->temp_photos,filename,pft->name,"temp_photos");
    fscanlimit2(file,&pft->temp,filename,pft->name,"temp");
    if(pft->type<0 || pft->type>=ntypes)
    {
      fprintf(stderr,"Invalid pft type=%d.\n",pft->type);
      break;
    }
    pft->mortality=nomortality;
    pft->init=noinit;
    pft->fire=nofire;
    pft->adjust=noadjust;
    if(scanfcn[pft->type](file,pft,filename))
      break;
  }
  pt_pclose(file);
  if(n<npft)
    fprintf(stderr,"Warning: number of PFTs truncated to %d.\n",n); 
  return n;
} /* of 'fscanpftpar' */
