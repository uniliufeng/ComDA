/***************************************************************************/
/**                                                                       **/
/**                    s  n  o  w  .  c                                   **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define tsnow 0.0
#define km 3.0
#define maxsnowpack 10000.0


/*    
 *    Function snow
 *
 *    Adjust daily precipitation by snowmelt and accumulation in snowpack
 *    Ref: Haxeltine & Prentice 1996
 *
 */

void snow(Real *snowpack, /* snowpack depth (mm) */
          Real *prec,     /* Precipitation (mm) */
          Real *snowmelt, /* snowmelt */
          Real temp       /* temperature */
         )
{
  if(temp<tsnow)
  {
    *snowpack+=*prec;
    if(*snowpack>maxsnowpack)
      *snowpack=maxsnowpack;
    *snowmelt=*prec=0.0;
  }
  else
  {
    *snowmelt=km*(temp-tsnow);
    if(*snowmelt>*snowpack) 
      *snowmelt=*snowpack;
    *snowpack-=*snowmelt;
  }
} /* of 'snow' */
