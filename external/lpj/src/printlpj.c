/***************************************************************************/
/**                                                                       **/
/**                      p  r  i  n  t  l  p  j  .  c                     **/
/**                                                                       **/
/**     print restart file of  LPJ                                        **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  20.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "grass.h"
#include "tree.h"

#define PRINTLPJ_VERSION "0.9.002"
#define NTYPES 2

int main(int argc,char **argv)
{
  int npft,nsoil,startgrid,ngridcell;
  Pftpar *pftpar;
  Soilpar *soilpar;
  Cell *grid;
  Config config;
  Fscanpftparfcn scanfcn[NTYPES]={fscanpft_grass,fscanpft_tree};
  printf("**** %s C Version %s (" __DATE__ ") ****\n",argv[0],PRINTLPJ_VERSION);
  if(fscanconfig(&config,&argc,&argv))
    fail("Error reading config\n");
  if(argc>0)
    startgrid=atoi(argv[0]);
  else
    startgrid=config.startgrid;
  if(argc>1)
    ngridcell=atoi(argv[1])-startgrid+1;
  else
    ngridcell=config.ngridcell;
  config.ngridcell=min(ngridcell,config.ngridcell-startgrid+config.startgrid); 
  if(startgrid>=config.startgrid)
    config.startgrid=startgrid;
  if((npft=fscanpftpar(&pftpar,config.pftpar_filename,scanfcn,NTYPES))==0)
    fail("Error reading '%s'",config.pftpar_filename);
  if((nsoil=fscansoilpar(&soilpar,config.soilpar_filename))==0)
    fail("Error reading '%s',",config.soilpar_filename);
  config.restart_filename=config.write_restart_filename;
  if((grid=newgrid(&config,pftpar,npft,soilpar,nsoil))==NULL)
    fail("Error initializing grid");
  printf("Year: %d\n",config.lastyear);
  printcell(grid,config.ngridcell);
  return EXIT_SUCCESS;
} /* of 'main' */
