/***************************************************************************/
/**                                                                       **/
/**                      g  e  t  c  l  i  m  a  t  e  .  c               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  30.05.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool getclimate(Climate *climate,int year)
{
  short *vec;
  int i;
  vec=newvec(short,climate->n);
  year-=climate->firstyear;
  fseek(climate->file_temp,year*climate->size+climate->offset,SEEK_SET);
  if(fread(vec,sizeof(short),climate->n,climate->file_temp)!=climate->n)
  {
    fprintf(stderr,"Error reading temperature of year %d in 'readclimate'.\n",year);
    free(vec);
    return TRUE;
  }
  if(climate->swap)
    for(i=0;i<climate->n;i++)
      climate->temp[i]=swapshort(vec[i])*0.1;
  else
    for(i=0;i<climate->n;i++)
      climate->temp[i]=vec[i]*0.1;
  fseek(climate->file_prec,year*climate->size+climate->offset,SEEK_SET);
  if(fread(vec,sizeof(short),climate->n,climate->file_prec)!=climate->n)
  {
    fprintf(stderr,"Error reading precipitation of year %d in 'readclimate'.\n",year);
    free(vec);
    return TRUE;
  }
  if(climate->swap)
    for(i=0;i<climate->n;i++)
      climate->prec[i]=swapshort(vec[i]);
  else
    for(i=0;i<climate->n;i++)
      climate->prec[i]=vec[i];
  fseek(climate->file_cloud,year*climate->size+climate->offset,SEEK_SET);
  if(fread(vec,sizeof(short),climate->n,climate->file_cloud)!=climate->n)
  {
    fprintf(stderr,"Error reading cloudiness of year %d in 'readclimate'.\n",year);
    free(vec);
    return TRUE;
  }
  if(climate->swap)
    for(i=0;i<climate->n;i++)
      climate->cloud[i]=100-swapshort(vec[i]);
  else
    for(i=0;i<climate->n;i++)
      climate->cloud[i]=100-vec[i];
  if(climate->file_wet!=NULL)
  {
    fseek(climate->file_wet,year*climate->size+climate->offset,SEEK_SET);
    if(fread(vec,sizeof(short),climate->n,climate->file_wet)!=climate->n)
    {
      fprintf(stderr,"Error reading wetdays of year %d in 'readclimate'.\n",year);
      free(vec);
      return TRUE;
    }
    if(climate->swap)
      for(i=0;i<climate->n;i++)
        climate->wet[i]=swapshort(vec[i]);
    else  
      for(i=0;i<climate->n;i++)
        climate->wet[i]=vec[i];
  }
  free(vec);
  return FALSE;
} /* of 'getclimate' */
