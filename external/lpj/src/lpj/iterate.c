/***************************************************************************/
/**                                                                       **/
/**                      i  t  e  r  a  t  e  .  c                        **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  30.05.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define nspinyear 30 /* cycle length during spinup (yr) */

int iterate(FILE **output,
            Cell *grid,
            Climate *climate,
            const Pftpar *pftpar,
            int npft,
            int ntypes,
            Config config)
{
  Real co2;
  int year;

  for(year=config.firstyear-config.nspinup;year<=config.lastyear;year++)
  {
    printf("Year: %d\n",year);
    if(iswriterestart(config) && year==config.lastyear)
     {
       printf("Write restart file '%s'.\n",config.write_restart_filename);
       fwriterestart(config.write_restart_filename,
		    grid,config.startgrid,config.ngridcell,
		    npft,year); /* write restart file */
     }
    co2=getco2(climate,year);
    if(getclimate(climate,(year<config.firstyear) ? (year-config.firstyear+config.nspinup) 
		  % nspinyear + config.firstyear : year))
      {
        fprintf(stderr,"Simulation stopped.\n");
        return year;
      }
    iterateyear(output,grid,climate,co2,pftpar,npft,ntypes,config,year); 
  } /* of 'for(year=...)' */
  return year;
} /* of 'iterate' */
