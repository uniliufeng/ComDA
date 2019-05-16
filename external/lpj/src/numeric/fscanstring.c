/***************************************************************************/
/**                                                                       **/
/**                  f  s  c  a  n  s  t  r  i  n  g  .  c                **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"

Bool fscanstring(FILE *file,char *s)
{
  char c;
  while(fscanf(file,"%c",&c)==1 && c!='\"');
  if(c!='\"')
    return TRUE; 
  for(; fscanf(file,"%c",s)==1 && (*s!='\"');s++);
  if(*s!='\"')
    return TRUE;
  *s='\0';
  return FALSE;
} /* of 'fscanstring' */
