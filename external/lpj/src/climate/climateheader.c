/***************************************************************************/
/**                                                                       **/
/**               c  l  i  m  a  t  e  h  e  a  d  e  r  .  c             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 19.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define LPJ_CLIMATE_HEADER "LPJCLIMATE"
#define LPJ_CLIMATE_VERSION 1

Bool fwriteclimateheader(FILE *file,Climateheader header)
{
  if(fwriteheader(file,LPJ_CLIMATE_HEADER,LPJ_CLIMATE_VERSION))
    return TRUE;
  if(fwrite(&header,sizeof(Climateheader),1,file)!=1)
    return TRUE;
  return FALSE;
} /* of 'fwriteclimateheader' */

int climateheadersize(void)
{
  return sizeof(int)+sizeof(Climateheader)+strlen(LPJ_CLIMATE_HEADER);
} /* of 'climateheadersize' */

Bool freadclimateheader(FILE *file,Climateheader *header,Bool *swap)
{
  if(freadheader(file,swap,LPJ_CLIMATE_HEADER,LPJ_CLIMATE_VERSION))
    return TRUE;
  if(fread(header,sizeof(Climateheader),1,file)!=1)
    return TRUE;
  if(*swap)
  {
    header->order=swapint(header->order);
    header->firstyear=swapint(header->firstyear);
    header->nyear=swapint(header->nyear);
    header->ncell=swapint(header->ncell);
  }
  return FALSE;
} /* of 'freadclimateheader' */
