/***************************************************************************/
/**                                                                       **/
/**                      u  n  i  t  s  .  h                              **/
/**                                                                       **/
/**     Header for unit conversions                                       **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 11.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef UNITS_H /* Already included? */
#define UNITS_H

/* Definition of macros */

#define deg2rad(deg) ((deg)*M_PI*0.00555555555555) /* Convert degree -> radian */
#define rad2deg(rad) ((rad)*180*M_1_PI) /* Convert radian -> degree */
#define degCtoK(deg) ((deg)+273.15)  /* Convert deg C --> Kelvin */
#define ppm2bar(ppm) ((ppm)*1e-6)      /* Convert ppmv --> bar */
#define ppm2Pa(ppm) ((ppm)*1e-1)      /* Convert ppmv --> Pa */
#define hour2sec(hour) ((hour)*3600)      /* Convert hour --> sec */
#define day2hour(day) ((day)*24)          /* Convert day --> hour */
#define hour2day(hour) ((hour)*0.041666666666666) /* Convert hour --> day */

#endif
