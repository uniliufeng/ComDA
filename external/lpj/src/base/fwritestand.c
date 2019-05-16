/***************************************************************************/
/**                                                                       **/
/**                 f  w  r  i  t  e  s  t  a  n  d  .  c                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.06.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

Bool fwritestand(FILE *file,const Stand *stand,Bool full)
{
  if(fwritepftlist(file,stand->pftlist,full)!=getnpft(stand->pftlist))
    return TRUE;
  if(fwritesoil(file,&stand->soil,full))
    return TRUE;
  return FALSE;
} /* of *fwritestand' */
