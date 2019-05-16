/***************************************************************************/
/**                                                                       **/
/**               i  n  i  t  c  l  i  m  a  t  e  .  c                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

static FILE *openfile(Climateheader *header,Bool *swap,const char *filename,
                      Config config)
{
  FILE *file;
  if((file=fopen(filename,"r"))==NULL)
  {
    printfopenerr("initclimate",filename);
    return NULL;
  }
  if(freadclimateheader(file,header,swap))
  {
    fprintf(stderr,"Invalid header in '%s'.\n",filename);
    fclose(file);
    return NULL;
  }
  if(header->firstyear>config.firstyear)
    fprintf(stderr,"Warning: first year in '%s'=%d greater than  %d.\n",
            filename,header->firstyear,config.firstyear);
  if(config.totalgridcell!=header->ncell)
    fprintf(stderr,"Warning: number of gridcells in '%s' distinct from %d.\n",
            filename,config.totalgridcell);
  return file;
} /* of 'openfile' */

/*
 *     Function initclimate
 *
 *     Initializes climate datatype 
 *
 */

Climate *initclimate(Config config)
{
  Climateheader header;
  FILE *file_co2;
  int i,yr;
  double val;
  Climate *climate;
  climate=new(Climate);
  if((climate->file_temp=openfile(&header,&climate->swap,config.temp_filename,
                                  config))==NULL)
  {
    free(climate);
    return NULL;
  }
  climate->firstyear=header.firstyear;
  if((climate->file_prec=openfile(&header,&climate->swap,config.prec_filename,
                                  config))==NULL)
  {
    fclose(climate->file_temp);
    free(climate);
    return NULL;
  }
  if((climate->file_cloud=openfile(&header,&climate->swap,
                                   config.cloud_filename,config))==NULL)
  {
    fclose(climate->file_temp);
    fclose(climate->file_prec);
    free(climate);
    return NULL;
  }
  if(!israndomprec(config))
    climate->file_wet=NULL;
  else 
  {
    if((climate->file_wet=openfile(&header,&climate->swap,config.wet_filename,
                                   config))==NULL)
    {
      fclose(climate->file_temp);
      fclose(climate->file_prec);
      fclose(climate->file_cloud);
      free(climate);
      return NULL;
    }
  }            
  if((file_co2=fopen(config.co2_filename,"r"))==NULL)
  {
    printfopenerr("initclimate",config.co2_filename);
    free(climate);
    return NULL;
  }
  climate->co2=newvec(Real,config.lastyear-config.firstyear+1+
                      config.firstyear-climate->firstyear);
  for(i=0;
      i<config.lastyear-config.firstyear+1+config.firstyear-climate->firstyear;
      i++)
    if(fscanf(file_co2,"%d %lg",&yr,&val)!=2)
    {
      fprintf(stderr,"Error reading co2 in 'readclimate'.\n");
      free(climate->co2);
      free(climate);
      return NULL;
    }
    else
    {
      if(i==0 && config.firstyear<yr)
        fprintf(stderr,"Warning: first year in '%s'=%d greater than %d.\n",
                config.co2_filename,yr,config.firstyear);
      climate->co2[i]=val;
    }
  fclose(file_co2);
//  if(header.order==CELLYEAR)
//  {
//    climate->size=config.totalgridcell*NMONTH*sizeof(short);
//    climate->n=NMONTH*config.ngridcell;
//    climate->offset=config.startgrid*NMONTH*sizeof(short)+climateheadersize();
//  }
//  else 
//  {
    climate->size=NMONTH*sizeof(short);
    climate->n=NMONTH*config.ngridcell; /*climate->n=header.nyear*NMONTH*config.ngridcell;*/
    climate->offset=config.startgrid*NMONTH*header.nyear*sizeof(short)+
                    climateheadersize();
//  }
  printf("climate->size:%d cliamte->n:%d, climate->offset%d \n",climate->size,climate->n,climate->offset);
  if((climate->cloud=newvec(Real,climate->n))==NULL)
  {
    printallocerr("initclimate","cloud");
    free(climate);
    return NULL;
  }
  if((climate->prec=newvec(Real,climate->n))==NULL)
  {
    printallocerr("initclimate","prec");
    free(climate->cloud);
    free(climate);
    return NULL;
  }
  if((climate->temp=newvec(Real,climate->n))==NULL)
  {
    printallocerr("initclimate","temp");
    free(climate->cloud);
    free(climate->prec);
    free(climate);
    return NULL;
  }
  if(israndomprec(config))
  {
    if((climate->wet=newvec(Real,climate->n))==NULL)
    {
      printallocerr("initclimate","wet");
      free(climate->co2);
      free(climate->cloud);
      free(climate->prec);
      free(climate->temp);
      free(climate);
      return NULL;
    }
  }
  return climate;
} /* of 'initclimate' */
