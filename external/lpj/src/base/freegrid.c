/***************************************************************************/
/**                                                                       **/
/**                    f  r  e  e  g  r  i  d  .  c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 22.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void freegrid(Cell *grid,int n)
{
  int cell;
  Stand *stand;
  for(cell=0;cell<n;cell++)
  {
    while(!isempty(grid[cell].standlist))
    {
       stand=getstand(grid[cell].standlist,0);
       freepftlist(stand->pftlist);
       free(stand);
       dellistitem(grid[cell].standlist,0);
    }
    freeclimbuf(&grid[cell].climbuf);
    free(grid[cell].gdd);
  }
  free(grid);
} /* of 'freegrid' */
