/***************************************************************************/
/**                                                                       **/
/**                    d  a  t  e  .  h                                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 23.09.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef DATE_H /* Already included? */
#define DATE_H

#define NMONTH 12 /* number of months in a year */
#define NDAYYEAR 365 /* number of days in a year */

typedef Real MReal[NMONTH];

extern int ndaymonth[NMONTH]; /* number of days in each month */
extern MReal ndaymonth1; /* 1/(number of days in each month) */
extern int midday[NMONTH+1];
extern MReal diffday;

/* Declaration of functions */

extern Real interpolate(const MReal,int,int);

/* Definition of macros */

#define foreachmonth(month) for(month=0;month<NMONTH;month++)
#define foreachdayofmonth(day,month) for(day=0;day<ndaymonth[month];day++)
#define foreachdayofyear(day) for(day=0;day<NDAYYEAR;day++)

#endif
