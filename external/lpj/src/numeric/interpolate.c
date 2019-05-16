/***************************************************************************/
/**                                                                       **/
/**                  i  n  t  e  r  p  o  l  a  t  e  .  c                **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"
#include "date.h"

/*
 *  Function interpol
 *
 *  Daily interpolation of monthly data
 *
 */

Real interpolate(const MReal mval, /* monthly values to interpolate       */
                 int month,        /* month (0..11)                       */
                 int dm            /* day of month (1 ..ndaymonth[month]) */
                )                  /* returns interpolated value          */
{
  int nextmonth,i;
  if(dm>ndaymonth[month]/2)
  {
    dm-=ndaymonth[month]/2;
    nextmonth=(month<NMONTH-1) ? month+1 : 0;
  }
  else
  {
    nextmonth=month;
    if(month==0)
      month=NMONTH-1;
    else
      month--; 
    dm+=ndaymonth[month]/2; 
  }
#ifdef DEBUG2
  printf("mval[%d]=%g,mval[%d]=%g,dm=%d,diffday=%g,%d\n",month,mval[month],
         nextmonth,mval[nextmonth],dm,diffday[month],ndaymonth[month]);
  printf("res=%g\n",mval[month]+dm*(mval[nextmonth]-mval[month])*diffday[month]);
#endif

  return mval[month]+dm*(mval[nextmonth]-mval[month])*diffday[month];
} /* of 'interpolate' */
