/***************************************************************************/
/**                                                                       **/
/**                   f  r  e  e  c  l  i  m  a  t  e  .  c               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

void freeclimate(Climate *climate)
{
  fclose(climate->file_temp);
  fclose(climate->file_prec);
  fclose(climate->file_cloud);
  free(climate->co2);
  if(climate->file_wet!=NULL)
  {
    fclose(climate->file_wet);
    free(climate->wet);
  }
  free(climate);
} /* of 'freeclimate' */
