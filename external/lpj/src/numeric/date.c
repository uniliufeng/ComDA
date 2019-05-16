/***************************************************************************/
/**                                                                       **/
/**                     d  a  t  e  .  c                                  **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  21.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"
#include "date.h"

int ndaymonth[NMONTH]=
{
  31,28,31,30,31,30,31,31,30,31,30,31
};

MReal ndaymonth1=
{
  1/31.0,1/28.0,1/31.0,1/30.0,1/31.0,1/30.0,1/31.0,1/31.0,1/30.0,1/31.0,
  1/30.0,1/31.0
};

int midday[NMONTH+1]=
{
  15,43,74,104,135,165,196,227,257,288,318,349,380 
};

MReal diffday=
{
  1/28.0,1/31.0,1/30.0,1/31.0,1/30.0,1/31.0,1/31.0,1/30.0,1/31.0,1/30.0,
  1/31.0,1/31.0
};
