/*
  Changes for the VIC implementation are preceded by the comment * start
  vic_change * and followed by the comment * end vic_change * */

/* RCS Id String
 * $Id: mtclim42_vic.h,v 5.4 2004/05/12 01:32:06 tbohn Exp $
 */

/* 
mtclim42.h
constants typedefs, and function prototypes for MTCLIM 4.2

Peter Thornton
NTSG, School of Forestry
University of Montana
5/10/98 

(dim) stands for dimensionless values

Adapted for inclusion in VIC-code:
Bart Nijssen
Sat Aug 21 16:58:43 1999
Last Changed: Fri Apr 25 11:37:15 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
*/



#define TDAYCOEF 0.45     /* daylight air temperature coefficient (dim) */

#define SECPERRAD 13750.9871     /* seconds per radian of hour angle */
#define RADPERDAY 0.017214       /* radians of Earth orbit per julian day */
#define RADPERDEG 0.01745329     /* radians per degree */
#define MINDECL -0.4092797       /* minimum declination (radians) */
#define DAYSOFF 11.25            /* julian day offset of winter solstice */
/* start vic_change */
#define SRADDT 30.0             /* timestep for radiation routine (seconds) */
				/* Note:  Make sure that 3600 % SRADDT == 0 */
/* end vic_change */

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define R        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
/* start vic_change */
#ifndef PI
#define PI       3.14159265
#endif
/* end vic_change */


