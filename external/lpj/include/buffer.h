/***************************************************************************/
/**                                                                       **/
/**                       b  u  f  f  e  r  .  h                          **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 12.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef BUFFER_H
#define BUFFER_H

/* Definitions of datatypes */

typedef struct
{
  int index,n,size;
  Real avg;
  Real *data;
} Buffer;

/* Declaration of functions */

extern Buffer *newbuffer(int);
extern Real updatebuffer(Buffer *,Real); 
extern Bool fwritebuffer(FILE *,const Buffer *);
extern Buffer *freadbuffer(FILE *,Bool);
extern void freebuffer(Buffer *);

#endif
