/***************************************************************************/
/**                                                                       **/
/**                c  a  t  l  p  j   .  c                                **/
/**                                                                       **/
/**     LPJ utility programme. Concatenates Restart files to stdout       **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 31.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include "types.h"
#include "swap.h"
#include "header.h"

#define USAGE "Usage: %s [-o filename] lpjfile ...\n"

typedef struct
{
  FILE *file;
  Restartheader header;
} Item;

int compare(const Item *a,const Item *b)
{
  if(a->header.firstcell<a->header.firstcell)
    return -1;
  else if(a->header.firstcell==a->header.firstcell)
    return 0;
  else
    return 1;
} /* of 'compare' */

int main(int argc,char **argv)
{
  int i,len,header_offset,o_offset,offset,j,ncell,swap,count;
  Item *item;
  int ptr;
  void *data;
  Restartheader header;
  struct stat filestat;
  FILE *out;
  out=stdout;
  for(i=1;i<argc;i++)
    if(argv[i][0]=='-')
    {
      if(!strcmp(argv[i],"-o"))
      {
        out=fopen(argv[++i],"wb");
        if(out==NULL)
        {
          fprintf(stderr,"Error creating '%s': %s\n",argv[i],strerror(errno));
          return EXIT_FAILURE;
        }
      }
      else
      {
        fprintf(stderr,"Invalid option '%s'.\n"
                USAGE,argv[i],argv[0]);
        return EXIT_FAILURE;
      }
    }
    else
      break;
  item=(Item *)malloc(sizeof(Item)*(argc-i));
  ncell=count=0;
  for(;i<argc;i++)
  {
    item[count].file=fopen(argv[i],"rb");
    if(item[count].file==NULL)
    {
      fprintf(stderr,"Error opening '%s': %s\n",argv[i],strerror(errno));
      continue;
    }
    if(freadheader(item[count].file,&swap,RESTART_HEADER,RESTART_VERSION))
    {
      fprintf(stderr,"Error reading header of file '%s'\n",argv[i]);
      continue;
    }
    if(swap)
    {
      fprintf(stderr,"Unsupported byte order in file '%s'\n",argv[i]);
      continue;
    }
    fread(&item[count].header,sizeof(Restartheader),1,item[count].file);
    ncell+=item[count].header.ncell;
    count++;
  }
  qsort(item,count,sizeof(Item),(int(*)(const void *,const void *))compare);
  header.year=item[0].header.year;
  for(i=1;i<count;i++)
  {
    if(item[i-1].header.firstcell+item[i-1].header.ncell!=item[i].header.firstcell) 
      fprintf(stderr,"Warning: first cell in %d=%d not last cell+1=%d\n",i,
              item[i].header.firstcell,
              item[i-1].header.firstcell+item[i-1].header.ncell);
    if(item[i].header.year!=header.year)
      fprintf(stderr,"Warning: year=%d in %d differs from %d\n",
              item[i].header.year,i,header.year);
  }
  header.firstcell=item[0].header.firstcell;
  header.year=item[0].header.year;
  header.ncell=ncell;
  fwriteheader(out,RESTART_HEADER,RESTART_VERSION);
  fwrite(&header,sizeof(header),1,out);
  header_offset=strlen(RESTART_HEADER)+sizeof(int)+sizeof(Restartheader);
  o_offset=header_offset;
  offset=0;
  for(i=0;i<count;i++)
  { 
    fstat(fileno(item[i].file),&filestat);
    len=filestat.st_size-header_offset-item[i].header.ncell*sizeof(int);
    fseek(out,o_offset,SEEK_SET);
    o_offset+=item[i].header.ncell*sizeof(int);
    for(j=0;j<item[i].header.ncell;j++)
    {
      fread(&ptr,sizeof(ptr),1,item[i].file);
      ptr+=(ncell-item[i].header.ncell)*sizeof(int)+offset;
      fwrite(&ptr,sizeof(ptr),1,out);
    }
    fseek(out,offset+header_offset+ncell*sizeof(int),SEEK_SET);
    offset+=len;
    data=malloc(len);
    fread(data,len,1,item[i].file);
    fwrite(data,len,1,out);
    free(data);
    fclose(item[i].file);
  }
  return EXIT_SUCCESS;
} /* of 'main' */
