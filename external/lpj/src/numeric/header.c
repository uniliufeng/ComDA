/***************************************************************************/
/**                                                                       **/
/**                     h  e  a  d  e  r  .  c                            **/
/**                                                                       **/
/**     Reading/Writing file header for  LPJ related files. Detects       **/
/**     whether byte order has to be changed                              **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  29.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "types.h"
#include "header.h"

Bool freadheader(FILE *file,Bool *swap,const char* header,int version)
{
  char *buffer;
  int file_version;
  buffer=(char *)malloc(strlen(header)+1);
  if(fread(buffer,strlen(header),1,file)!=1)
  {
    free(buffer);
    return TRUE;
  }
  buffer[strlen(header)]='\0';
  if(strcmp(buffer,header))
  {
    free(buffer);
    return TRUE;
  }
  free(buffer);
  if(fread(&file_version,sizeof(file_version),1,file)!=1)
    return TRUE;
  if((file_version & 0xff)==0)
  {
    *swap=TRUE;
    file_version=swapint(file_version);
  }
  else
    *swap=FALSE;
  return (file_version!=version);
} /* of 'freadheader' */

Bool fwriteheader(FILE *file,const char *header,int version)
{
  if(fwrite(header,strlen(header),1,file)!=1)
    return TRUE;
  return (fwrite(&version,sizeof(version),1,file)!=1);
} /* of 'fwriteheader' */
