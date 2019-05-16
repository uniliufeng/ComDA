/***************************************************************************/
/**                                                                       **/
/**                      c  r  u  2  c  l  m  .  c                        **/
/**                                                                       **/
/**                                                                       **/
/**     Converts CRU data files into LPJ climate data files               **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 19.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "types.h"
#include "swap.h"
#include "date.h"
#include "buffer.h"
#include "coord.h"
#include "config.h"
#include "climate.h"

#define CRU2CLM_VERSION "1.0.001"
#define NYEAR 103  /* default value for number of years */
#define FIRSTYEAR 1901 /* default value for first year */
#define FIRSTCELL 0
#define NCELL 67420
#define USAGE    "Usage: cru2clm [-h] [-firstyear first] [-lastyear last]\n"\
                 "       [-nyear n] [-firstcell] [-ncell] [-nbands] [-rev]\n"\ 
                 "       [-yearcell] crufile clmfile\n"

typedef short Data[NMONTH];

int main(int argc,char **argv)
{
  FILE *file;
  Data *data;
  Climateheader header;
  long int k;
  int i,j,n;
  Bool revorder;
  struct stat filestat;
  header.nyear=NYEAR;
  header.firstyear=FIRSTYEAR;
  header.order=CELLYEAR;
  header.ncell=NCELL;
  revorder=FALSE;
  for(i=1;i<argc;i++)
    if(argv[i][0]=='-')
    {
      if(!strcmp(argv[i],"-h"))
      {
        printf("cru2clm " CRU2CLM_VERSION " (" __DATE__ ") - convert cru data files to\n"
               "       clm data files for lpj C version\n");
        printf(USAGE
               "-h               print this help text\n" 
               "-firstyear first first year in cru file (default %d)\n"
               "-lastyear last   last year in cru file\n"
               "-nyear n         number of years in cru file (default %d)\n"
               "-rev             change byte order in cru file\n"
               "-yearcell        does not revert order in cru file\n"
               "crufile          filename of cru data file"
               "clmfile          filename of clm data file",
               FIRSTYEAR,NYEAR);
        return EXIT_SUCCESS;
      }
      else if(!strcmp(argv[i],"-nyear"))
        header.nyear=atoi(argv[++i]);
      else if(!strcmp(argv[i],"-firstyear"))
        header.firstyear=atoi(argv[++i]);
      else if(!strcmp(argv[i],"-lastyear"))
        header.nyear=atoi(argv[++i])-header.firstyear+1;
      else if(!strcmp(argv[i],"-rev"))
        revorder=TRUE;
      else if(!strcmp(argv[i],"-yearcell"))
        header.order=YEARCELL;
      else
      {
        fprintf(stderr,"Error: invalid option '%s'.\n",argv[i]);
        fprintf(stderr,USAGE);
        return EXIT_FAILURE;
      }
    }
    else 
      break;
  if(argc<i+2)
  {
    fprintf(stderr,"Filenames missing.\n");
    fprintf(stderr,USAGE);
    return EXIT_FAILURE;
  }
  file=fopen(argv[i],"r");
  if(file==NULL)
  {
    fprintf(stderr,"Error opening '%s': %s\n",argv[i],strerror(errno));
    return EXIT_FAILURE;
  }
  fstat(fileno(file),&filestat);
  n=filestat.st_size/header.nyear/sizeof(Data);
  printf("Number of cells: %d\n",n);
  data=(Data *)malloc(n*sizeof(Data)*header.nyear);
  fread(data,sizeof(Data),n*header.nyear,file);
  if(revorder)
    for(k=0;k<n*header.nyear;k++)
      for(j=0;j<NMONTH;j++)
        data[k][j]=swapshort(data[k][j]); 
  fclose(file);
  file=fopen(argv[i+1],"w");
  if(file==NULL)
  {
    fprintf(stderr,"Error creating '%s': %s\n",argv[i+1],strerror(errno));
    return EXIT_FAILURE;
  }
  header.ncell=n;
  fwriteclimateheader(file,header);
  if(header.order==CELLYEAR)
  {
    for(i=0;i<header.nyear;i++)
      for(j=0;j<n;j++)
        fwrite(data+j*header.nyear+i,sizeof(Data),1,file);
  }
  else
    fwrite(data,sizeof(Data),header.nyear*n,file);
  free(data);
  fclose(file);
  return EXIT_SUCCESS;
} /* of 'main' */
