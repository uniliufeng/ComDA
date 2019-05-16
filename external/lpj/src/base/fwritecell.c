/***************************************************************************/
/**                                                                       **/
/**             f  w  r  i  t  e  c  e  l  l  .  c                        **/
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

int fwritecell(FILE *file,int index[],const Cell grid[],int ncell,int npft,
               Bool full)
{
  int cell;
  for(cell=0;cell<ncell;cell++)
  {
    if(index!=NULL)
      index[cell]=ftell(file);
    fwrite(&grid[cell].skip,sizeof(Bool),1,file);
    if(!grid[cell].skip)
    {
      if(fwritestandlist(file,grid[cell].standlist,full)!=
         grid[cell].standlist->n)
        break;
      if(full)
      {
        if(fwriteclimbuf(file,&grid[cell].climbuf))
          break;
        if(fwrite(grid[cell].gdd,sizeof(Real),npft,file)!=npft)
          break;
      }
    }
  }
  return cell;
} /* of 'fwritecell' */
