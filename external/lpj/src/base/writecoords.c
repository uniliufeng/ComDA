/***************************************************************************/
/**                                                                       **/
/**                   w  r  i  t  e  c  o  o  r  d  s  .  c               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 04.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

int writecoords(FILE *file,const Cell grid[],int ncell)
{
  int cell;
  for(cell=0;cell<ncell;cell++)
    if(!grid[cell].skip)
    {
      if(writecoord(file,grid[cell].coord))
        break;
    }
  return cell;
} /* of 'writecoords' */
