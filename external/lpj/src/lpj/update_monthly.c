/***************************************************************************/
/**                                                                       **/
/**              u  p  d  a  t  e  _  m  o  n  t  h  l  y  .  c           **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"


Real update_monthly( FILE **output,
		     Cell *cell,
		     Real mtemp,
		     int month,
		     int year,
		     int npft,
		     Config config)
{
  int s,l;
  Real nep;
  Stand *stand;
  monthly_climbuf(&cell->climbuf,mtemp);
  foreachstand(stand,s,cell->standlist)
    getlag(&stand->soil,month);
  for(l=0;l<NSOILLAYER;l++)
    cell->output.mswc[l]*=ndaymonth1[month];
  if(year>=config.firstyear) 
    fwriteoutput_monthly(output,cell,config.n_out,npft);
  nep=cell->output.mnpp-cell->output.mrh;
  initoutput_monthly(&(cell->output));
  return nep;
} /* of 'monthly_update' */  
