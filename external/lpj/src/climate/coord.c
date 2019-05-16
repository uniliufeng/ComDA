/***************************************************************************/
/**                                                                       **/
/**                     c  o  o  r  d  .  c                               **/
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

#define LPJGRID_HEADER "LPJGRID"
#define LPJGRID_VERSION 1

typedef struct
{
  short lon,lat;
} Icoord;

Coordfile *opencoord(const char *filename /* filename of coord file */
                    )           
{
  Coordfile *coordfile;
  coordfile=new(Coordfile);
  coordfile->file=fopen(filename,"r");
  if(coordfile->file==NULL)
  {
    free(coordfile);
    return NULL;
  }
  if(freadheader(coordfile->file,&coordfile->swap,LPJGRID_HEADER,LPJGRID_VERSION))
  {
    fclose(coordfile->file);
    free(coordfile);
    return NULL;
  }
  fread(&coordfile->n,sizeof(int),1,coordfile->file);
  if(coordfile->swap)
    coordfile->n=swapint(coordfile->n);
  return coordfile;  
} /* of 'opencoord' */

void closecoord(Coordfile *coordfile)
{
  fclose(coordfile->file);
  free(coordfile);
} /* of 'closecoord' */

Bool readcoord(Coordfile *coordfile,
               Coord *coord,         /* cell coordinate read from file */
               Coord resol           /* resolution (deg) */
              )                      /* returns FALSE for successful read */
{
  Icoord icoord;
  if(fread(&icoord,sizeof(icoord),1,coordfile->file)!=1)
    return TRUE;
  if(coordfile->swap)
  { 
    coord->lat=swapshort(icoord.lat)*0.01;
    coord->lon=swapshort(icoord.lon)*0.01;
  }
  else
  { 
    coord->lat=icoord.lat*0.01;
    coord->lon=icoord.lon*0.01;
  }
  coord->area=cellarea(*coord,resol);
  return FALSE;
} /* of 'getcoord' */

Bool writecoord(FILE *file,  /* file pointer */
                Coord coord  /* cell coordinat written to file */
               )             /* returns FALSE for successful write */
{
  Icoord icoord;
  icoord.lat=(int)(coord.lat*100);
  icoord.lon=(int)(coord.lon*100);
  return fwrite(&icoord,sizeof(icoord),1,file)!=1;
} /* of 'writecoord' */

int seekcoord(Coordfile *coordfile, /* file pointer */
              int pos               /*  position in file */
             )                      /* returns return code of fseek */
{
  return fseek(coordfile->file,
               pos*sizeof(Icoord)+strlen(LPJGRID_HEADER)+2*sizeof(int),
               SEEK_SET);
} /* of 'seekcoord' */

Real cellarea(Coord coord, /* cell coordnate */
              Coord resol  /* Resolution (deg) */
             )             /* returns area of cell (m^2) */
{
  return (111e3*resol.lat)*(111e3*resol.lon)*cos(deg2rad(coord.lat));
} /* of 'cellarea' */
