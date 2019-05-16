/***************************************************************************/
/**                                                                       **/
/**                      i  t  e  r  a  t  e  y  e  a  r  .  c            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  22.11.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define nspinyear 30 /* cycle length during spinup (yr) */

void iterateyear(FILE **output,
		 Cell *grid,
		 Climate *climate,
		 Real co2,
		 const Pftpar *pftpar,
		 int npft,
		 int ntypes,
		 Config config,
		 int year)
{
  Real temp,prec,sun;

  int month,dayofmonth,day;
  int cell;
  Real *nep,balanceC,delta_totc;
  
  nep=newvec(Real,config.ngridcell);
  
  for(cell=0;cell<config.ngridcell;cell++){  
    if(!(grid+cell)->skip)
    {
      initoutput_annual(&((grid+cell)->output));
      init_annual(grid+cell,npft);
      nep[cell]=0.0;
    }
  }

  day=1;
  foreachmonth(month){
    for(cell=0;cell<config.ngridcell;cell++)
      if(!(grid+cell)->skip){
	if(israndomprec(config))
	  prdaily(&((grid+cell)->climbuf),ndaymonth[month],
		  (getcellprec(climate,cell))[month],
		  (getcellwet(climate,cell))[month]);

#ifdef DEBUG
      printf("temp = %.2f prec = %.2f sun = %.2f wet = %.2f\n",
	     (getcelltemp(climate,cell))[month],
	     (getcellprec(climate,cell))[month],
	     (getcellsun(climate,cell))[month],
	     (israndomprec(config)) ? (getcellwet(climate,cell))[month] : 0);
#endif 
}
    foreachdayofmonth(dayofmonth,month){
      for(cell=0;cell<config.ngridcell;cell++){
	if(!(grid+cell)->skip){
	
	  temp=interpolate(getcelltemp(climate,cell),month,dayofmonth);
	  sun=interpolate(getcellsun(climate,cell),month,dayofmonth);
	  prec=(israndomprec(config)) ? getprec(*(grid+cell),dayofmonth) :
	    interpolate(getcellprec(climate,cell),month,dayofmonth)*
	    ndaymonth1[month];
	  prec=prec>0.000001 ? prec : 0.0;
	  updategdd(grid[cell].gdd,pftpar,npft,temp);
	  update_daily(grid+cell,co2,temp,prec,sun,day);
	  (grid+cell)->aprec+=prec;
	}
      }
      day++;
    }/* of 'foreachdayofmonth */

    for(cell=0;cell<config.ngridcell;cell++)
      if(!(grid+cell)->skip)
	nep[cell]+=update_monthly(output,grid+cell,(getcelltemp(climate,cell))[month],month,year,npft,config);
  } /* of 'foreachmonth */

  for(cell=0;cell<config.ngridcell;cell++){
    if(!(grid+cell)->skip){
      delta_totc=update_annual(output,grid+cell,pftpar,npft,ntypes,year,config);
      
      balanceC=nep[cell]-grid[cell].output.firec+grid[cell].output.flux_estab-delta_totc;
      if (year>config.firstyear)
       if (balanceC>2 || balanceC<-2)
         fail("BALANCE-error %.2f\n",balanceC);
      
#ifdef DEBUG
      printf("year=%d\n",year);
      printf("cell=%d\n",cell);
      printcell(grid+cell,1);
#endif
      if(config.nspinup>soil_equil_year && 
	 year==config.firstyear-config.nspinup+soil_equil_year)
	equilsom(grid+cell);
    }
  }
  free(nep);
} /* of 'iterateyear' */
