/***************************************************************************/
/**                                                                       **/
/**                     b  u  f  f  e  r  .  c                            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "types.h"
#include "swap.h"
#include "errmsg.h"
#include "buffer.h"

Buffer *newbuffer(int size)
{
  Buffer *buffer;
  buffer=new(Buffer);
  buffer->size=size;
  buffer->n=buffer->index=0;
  if((buffer->data=newvec(Real,size))==NULL)
  {
    printallocerr("initbuffer","data");
    free(buffer);
    return NULL;
  }
  buffer->avg=0;
  return buffer;
} /* of 'newbuffer' */

Real updatebuffer(Buffer *buffer,Real val)
{
  if(buffer->n<buffer->size)
  {
    buffer->data[buffer->n]=val;
    buffer->n++;
    buffer->avg+=val;
  }
  else
  {
    buffer->avg-=buffer->data[buffer->index];
    buffer->data[buffer->index]=val;
    buffer->index=(buffer->index+1)% buffer->size;
    buffer->avg+=val;
  }
  return buffer->avg/buffer->n;
} /* of 'updatebuffer' */

Bool fwritebuffer(FILE *file,const Buffer *buffer)
{
  fwrite(&buffer->size,sizeof(int),1,file);
  fwrite(&buffer->n,sizeof(int),1,file);
  fwrite(&buffer->index,sizeof(int),1,file);
  fwrite(&buffer->avg,sizeof(Real),1,file);
  return (fwrite(buffer->data,sizeof(Real),buffer->n,file)!=buffer->n);
} /* of 'fwritebuffer' */

Buffer *freadbuffer(FILE *file,Bool swap)
{ 
  int size;
  Buffer *buffer;
  if(freadint1(&size,swap,file)!=1)
    return NULL;
  buffer=new(Buffer);
  buffer->size=size;
  buffer->data=newvec(Real,buffer->size);
  freadint1(&buffer->n,swap,file);
  freadint1(&buffer->index,swap,file);
  freadreal1(&buffer->avg,swap,file);
  if(freadreal(buffer->data,buffer->n,swap,file)!=buffer->n)
  {
    /* read error occured */
    free(buffer->data);
    free(buffer);
    return NULL;
  }
  return buffer;
} /* of 'freadbuffer' */

void freebuffer(Buffer *buffer)
{
  free(buffer->data);
  free(buffer);
} /* of 'freebuffer' */
