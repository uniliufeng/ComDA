/***************************************************************************/
/**                                                                       **/
/**               f  p  r  i  n  t  c  e  l  l  .  c                      **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void fprintcell(FILE *file,const Cell grid[],int ncell)
{
  int cell;
  for(cell=0;cell<ncell;cell++)
  {
    fprintf(file,"Coord :");
    fprintcoord(file,grid[cell].coord);
    fprintf(file,"\n");
    if(grid[cell].skip)
      fprintf(file,"Invalid soil\n");
    else
      fprintstandlist(file,grid[cell].standlist);
  }
} /* of 'fprintcell' */ 
