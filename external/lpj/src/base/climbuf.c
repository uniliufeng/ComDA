/***************************************************************************/
/**                                                                       **/
/**                 c  l  i  m   b  u  f  .  c                            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 09.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define k (1.0/12.0)

void new_climbuf(Climbuf *climbuf)
{
  int d;
  climbuf->max=newbuffer(CLIMBUFSIZE);
  climbuf->min=newbuffer(CLIMBUFSIZE);
  climbuf->atemp_mean=0;
  climbuf->dval_prec[0]=0;
  for(d=0;d<NDAYS;d++)
    climbuf->temp[d]=0;
} /* of 'new_climbuf' */

void init_climbuf(Climbuf *climbuf)
{
  climbuf->temp_min=HUGE_VAL;
  climbuf->temp_max=-HUGE_VAL;
  climbuf->gdd5=0;
} /* of 'init_climbuf' */

void daily_climbuf(Climbuf *climbuf,Real temp)
{
  int d;
  updategdd5(climbuf,temp);
  for(d=1;d<NDAYS;d++)
    climbuf->temp[d-1]=climbuf->temp[d];
  climbuf->temp[NDAYS-1]=temp;
} /* of  'daily_climbuf' */

void monthly_climbuf(Climbuf *climbuf,Real mtemp)
{
  int d;
  Real temp_avg;
  temp_avg=0;
  for(d=0;d<NDAYS;d++)
    temp_avg+=climbuf->temp[d];
  temp_avg/=NDAYS;  
  if(climbuf->temp_min>mtemp)
    climbuf->temp_min=mtemp;
  if(climbuf->temp_max<mtemp)
    climbuf->temp_max=mtemp;
  climbuf->atemp_mean=(1-k)*climbuf->atemp_mean+k*mtemp;
} /* of 'monthly_climbuf' */
  
void annual_climbuf(Climbuf *climbuf)
{
  updatebuffer(climbuf->min,climbuf->temp_min);
  updatebuffer(climbuf->max,climbuf->temp_max);  
} /* of 'annual_climbuf' */

Bool fwriteclimbuf(FILE *file,const Climbuf *climbuf)
{
  fwrite(&climbuf->temp_max,sizeof(Real),1,file);
  fwrite(&climbuf->temp_min,sizeof(Real),1,file);
  fwrite(&climbuf->atemp_mean,sizeof(Real),1,file);
  fwrite(&climbuf->gdd5,sizeof(Real),1,file);
  fwrite(climbuf->temp,sizeof(Real),NDAYS,file);
  fwritebuffer(file,climbuf->min);
  return fwritebuffer(file,climbuf->max);
} /* of 'fwriteclimbuf' */

Bool freadclimbuf(FILE *file,Climbuf *climbuf,Bool swap)
{
  freadreal1(&climbuf->temp_max,swap,file);
  freadreal1(&climbuf->temp_min,swap,file);
  freadreal1(&climbuf->atemp_mean,swap,file);
  freadreal1(&climbuf->gdd5,swap,file);
  freadreal(climbuf->temp,NDAYS,swap,file);
  climbuf->min=freadbuffer(file,swap);
  climbuf->max=freadbuffer(file,swap);
  return (climbuf->max==NULL);
} /* of 'freadclimbuf' */

void freeclimbuf(Climbuf *climbuf)
{
  free(climbuf->max);
  free(climbuf->min);
} /* of 'freeclimbuf' */
