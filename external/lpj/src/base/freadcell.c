/***************************************************************************/
/**                                                                       **/
/**             f  r  e  a  d  c  e  l  l  .  c                           **/
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

#include "lpj.h"

Bool freadcell(FILE *file,Cell *cell,const Pftpar pftpar[],int npft,
               const Soilpar *soilpar,Bool swap)
{
  freadint(&cell->skip,1,swap,file);
  if(!cell->skip)
  {
    cell->standlist=freadstandlist(file,pftpar,soilpar,swap);
    if(cell==NULL)
      return TRUE;
    freadclimbuf(file,&cell->climbuf,swap);
    freadreal(cell->gdd,npft,swap,file);
  }
  return FALSE;
} /* of 'freadcell' */
