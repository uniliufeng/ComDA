/***************************************************************************/
/**                                                                       **/
/**                         l  p  j  .  c                                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  27.05.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#include <time.h>
#ifdef USE_UNAME
#include <sys/utsname.h>
#endif
#include "lpj.h"
#include "grass.h"
#include "tree.h"

#define NTYPES 2 /* number of PFT types: grass, tree */

int main(int argc,char **argv)
{
  FILE **output;
  char *progname;
  int npft,nsoil,year;
  Pftpar *pftpar;
  Soilpar *soilpar;
  Cell *grid;
  Climate *climate;
  time_t tstart,tend;
#ifdef USE_UNAME
  struct utsname osname;
#endif
  Fscanpftparfcn scanfcn[NTYPES]={fscanpft_grass,fscanpft_tree};
  Config config;
#ifdef USE_UNAME
  uname(&osname);
  printf("**** %s C Version %s (" __DATE__ ") %s %s.%s ****\n\n",argv[0],
         LPJ_VERSION,osname.sysname,osname.version,osname.release);
#else
  printf("**** %s C Version %s (" __DATE__ ") ****\n\n",argv[0],
         LPJ_VERSION);
#endif
  progname=argv[0];
  if(fscanconfig(&config,&argc,&argv))
    fail("Error reading config");
  if((npft=fscanpftpar(&pftpar,config.pftpar_filename,scanfcn,NTYPES))==0)
    fail("Error reading '%s'",config.pftpar_filename);
  if((nsoil=fscansoilpar(&soilpar,config.soilpar_filename))==0)
    fail("Error reading '%s'",config.soilpar_filename);
  if(isreadrestart(config))
    printf("Starting from restart file '%s'.\n",config.restart_filename);
  if((grid=newgrid(&config,pftpar,npft,soilpar,nsoil))==NULL)
    fail("Error initializing grid");
  if((climate=initclimate(config))==NULL) 
    fail("Error initializing climate");
  /* Initialize random seed */
#ifdef USE_RAND48
  srand48((config.seed==RANDOM_SEED) ? time(NULL) : config.seed);  
#else
  setseed((config.seed==RANDOM_SEED) ? time(NULL) : config.seed);  
#endif
  output=fopenoutput(config);
  printf("Simulation begins...\n");
  time(&tstart); /* Start timing */
  year=iterate(output,grid,climate,pftpar,npft,NTYPES,config);
  time(&tend); /* Stop timing */
  if(output[GRID]!=NULL)
    writecoords(output[GRID],grid,config.ngridcell);
  fcloseoutput(output,config.n_out);
  printf("Simulation ended.\n");
  /* free memory */
  freeclimate(climate);
  freegrid(grid,config.ngridcell);
  printf( (year>config.lastyear) ? "%s" : "%s successfully",progname);
  printf(" terminated, %d grid cells processed.\n"
         "Wall clock time: %d sec, %.2g sec/cell/year.\n",
         config.ngridcell,(int)(tend-tstart),
         (double)(tend-tstart)/config.ngridcell/(year-config.firstyear+
                                                 config.nspinup));
  return EXIT_SUCCESS;
} /* of 'main' */
