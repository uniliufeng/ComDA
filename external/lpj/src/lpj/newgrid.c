/***************************************************************************/
/**                                                                       **/
/**                       n  e  w  g  r  i  d  .  c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 09.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Cell *newgrid(Config *config, /* Pointer to LPJ configuration */
              const Pftpar pftpar[],
              int npft,       /* number of PFT's */
              const Soilpar soilpar[],
              int nsoil      /* number of soil types */
             )               /* returns allocated cell grid */
{
  Cell *grid;
  Stand *stand;
  Restartheader header;
  int i,offset;
  Bool swap;
  char soilcode;
  FILE *file_soil,*file_restart;
  Coordfile *file_coord;
  file_coord=opencoord(config->coord_filename);
  if(file_coord==NULL)
  {
    printfopenerr("newgrid",config->coord_filename);
    return NULL;
  }
  file_soil=fopen(config->soil_filename,"r");
  if(file_soil==NULL)
  {
    printfopenerr("newgrid",config->soil_filename);
    return NULL;
  }
  config->totalgridcell=numcoord(file_coord);
  if(config->ngridcell==ALL)
    config->ngridcell=config->totalgridcell;
  if(config->restart_filename!=NULL)
  {
    file_restart=fopen(config->restart_filename,"r");
    if(file_restart==NULL)
    {
      printfopenerr("newgrid",config->restart_filename);
      return NULL;
    }
    if(freadheader(file_restart,&swap,RESTART_HEADER,RESTART_VERSION))
    {
      fprintf(stderr,"Invalid header in '%s'.\n",config->restart_filename);
      fclose(file_restart);
      return NULL;
    }
    if(fread(&header,sizeof(header),1,file_restart)!=1)
    {
      fprintf(stderr,"Error reading '%s'.\n",config->restart_filename);
      fclose(file_restart);
      return NULL;
    }
    if(swap)
    {
      header.year=swapint(header.year);
      header.firstcell=swapint(header.firstcell);
      header.ncell=swapint(header.ncell);
    }
    if(header.year!=config->firstyear)
      fprintf(stderr,"Warning: year of restartfile=%d not equal start year=%d.\n",header.year,config->firstyear);
 
    offset=config->startgrid-header.firstcell;
    if(offset<0)
    {
      fprintf(stderr,"Warning: first grid cell not in restart file, set to %d.\n",header.firstcell);
      config->startgrid=header.firstcell;
      offset=0;
    }
    if(offset>0)
    {
      if(offset>=header.ncell)
      {
        fprintf(stderr,"Error: First grid cell not in restart file.\n");
        fclose(file_restart);
        return NULL;
      }
      fseek(file_restart,offset*sizeof(int),SEEK_CUR);
    }
    freadint1(&offset,swap,file_restart);
    fseek(file_restart,offset,SEEK_SET);
    if(header.ncell<config->ngridcell)
    {
      fprintf(stderr,"Warning: restart file too short, grid truncated to %d.\n",              header.ncell);
      config->ngridcell=header.ncell;
    }
           
  } 
  /* allocate grid */ 
  if((grid=newvec(Cell,config->ngridcell))==NULL)
  {
    printallocerr("newgrid","grid");
    return NULL;
  }
  seekcoord(file_coord,config->startgrid);
  fseek(file_soil,config->startgrid*sizeof(soilcode),SEEK_SET);
  for(i=0;i<config->ngridcell;i++)
  {
    if(readcoord(file_coord,&grid[i].coord,config->resolution))
    {
      fprintf(stderr,"Unexpected end of file in '%s', number of gridcells truncated to %d.\n",config->coord_filename,i);
      config->ngridcell=i;
      break;
    } 
    if(fread(&soilcode,sizeof(soilcode),1,file_soil)!=1)
    {
      fprintf(stderr,"Unexpected end of file in '%s', number of gridcells truncated to %d.\n",config->soil_filename,i);
      config->ngridcell=i;
      break;
    }
    grid[i].gdd=newgdd(npft);
    if(config->restart_filename==NULL)
    {
      grid[i].standlist=newlist(); 
      stand=new(Stand);
      addlistitem(grid[i].standlist,stand);
      stand->frac=1;
      stand->pftlist=newpftlist();
      if(soilcode<1 || soilcode>nsoil)
      {
        fprintf(stderr,"Invalid soilcode=%d, cell skipped\n",soilcode);
        grid[i].skip=TRUE;
      }
      else
      {
        initsoil(&stand->soil,soilpar+soilcode-1);
        grid[i].skip=FALSE;
      }
      new_climbuf(&grid[i].climbuf);
    }
    else if(freadcell(file_restart,grid+i,pftpar,npft,soilpar+soilcode-1,swap))
    {
      fprintf(stderr,"Unexpected end of file in '%s', number of gridcells truncated to %d.\n",config->restart_filename,i);
      config->ngridcell=i;
      break;
    }
  
  } /* of for(i=0;...) */
  if(config->restart_filename!=NULL)
    fclose(file_restart);
  closecoord(file_coord); 
  fclose(file_soil); 
  return grid;
} /* of 'newgrid' */
