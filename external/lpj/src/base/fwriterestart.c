/***************************************************************************/
/**                                                                       **/
/**                 f  w  r  i  t  e  r  e  s  t  a  r  t  .  c           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 29.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool fwriterestart(const char *name, /* filename                 */
                   const Cell *grid, /* cell array               */
                   int firstcell,    /* index of first cell      */
                   int ncell,        /* number of cells          */
                   int npft,         /* number of PFT's          */
                   int year          /* year                     */
                  )                  /* returns FALSE on success */ 
{
  FILE *file;
  int *index;
  int pos;
  Restartheader header;
  file=fopen(name,"w");
  if(file==NULL)
  {
    printfopenerr("fwriterestart",name);   
    return TRUE;
  }
  index=newvec(int,ncell);
  fwriteheader(file,RESTART_HEADER,RESTART_VERSION);
  header.year=year;
  header.firstcell=firstcell;
  header.ncell=ncell;
  fwrite(&header,sizeof(header),1,file);
  pos=ftell(file);
  fseek(file,sizeof(int)*ncell,SEEK_CUR);
  if(fwritecell(file,index,grid,ncell,npft,TRUE)!=ncell)
  {
    fprintf(stderr,"Error writing restart file '%s': %s\n",name,
            strerror(errno));
    free(index);
    fclose(file);
    return TRUE;
  }
  fseek(file,pos,SEEK_SET);
  fwrite(index,sizeof(int),ncell,file);
  fclose(file);
  free(index);
  return FALSE;
} /* of 'fwriterestart' */
