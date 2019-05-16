/***************************************************************************/
/**                                                                       **/
/**                         r  a  n  d  .  c                              **/
/**                                                                       **/
/**     Random number generator                                           **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 03.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include "types.h"
#include "numeric.h"

#define modulus 2147483647
#define multiplier 16807 
#define q 127773
#define r 2836

static long seed=12345678; /* seed for random number generator (see randfrac) */

void setseed(int init) 
{
  seed=init;
} /* of 'setseed' */

int getseed(void)
{
  return seed;
} /* of 'getseed' */

Real randfrac(void) 
{

 /*
  * DESCRIPTION
  * Returns a random floating-point number in the range 0-1.
  * Uses and updates the global variable 'seed' which may be initialised to any
  * positive integral value (the same initial value will result in the same 
  * sequence of returned values on subsequent calls to randfloat every time
  * the program is run)
  *
  * Reference: Park & Miller 1988 CACM 31: 1192
  *
  */

  seed=multiplier*(seed%q)-r*seed/q;
  if (!seed) 
    seed++; /* increment seed to 1 in unlikely event of 0 value */
  else if (seed<0) 
    seed+=modulus;
  return (Real)seed/(Real)modulus;
} /* of 'randfrac' */
