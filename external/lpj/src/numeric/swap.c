/***************************************************************************/
/**                                                                       **/
/**                     s  w  a  p  .  c                                  **/
/**                                                                       **/
/**     Converts big endian into littel endian                            **/
/**     Needed for reading data under Linux                               **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  14.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "swap.h"

typedef struct
{
  int lo,hi;
}Num;

static void swap(char *a,char *b)
{
  char h;
  h=*a;
  *a=*b;
  *b=h;
} /* of 'swap' */

short int swapshort(short int x)
{
  swap((char *)&x,(char *)(&x)+1);
  return x;
} /* of 'swapshort' */

int swapint(int x)
{
  swap((char *)&x,(char *)(&x)+3);
  swap((char *)&x+1,(char *)(&x)+2);
  return x;
} /* of 'swapint' */

static double swapdouble(Num num)
{
  double ret;
  Num x;
  x.hi=swapint(num.lo);
  x.lo=swapint(num.hi);
  memcpy(&ret,&x,sizeof(Num));
  return ret;
} /* of 'swapdouble' */

int freadreal(Real *data,int n,Bool swap,FILE *file)
{
  Num *num;
  int i,rc;
  if(swap)
  {
    num=(Num *)malloc(sizeof(Num)*n);
    rc=fread(num,sizeof(num),n,file); 
    for(i=0;i<rc;i++)
      data[i]=swapdouble(num[i]);
    free(num);
  }
  else
    rc=fread(data,sizeof(double),n,file);
  return rc;
} /* of 'freadreal' */

int freadint(int *data,int n,Bool swap,FILE *file)
{
  int i,rc;
  rc=fread(data,sizeof(int),n,file);
  if(swap)
    for(i=0;i<rc;i++)
      data[i]=swapint(data[i]);
  return rc;
} /* of 'freadint' */
