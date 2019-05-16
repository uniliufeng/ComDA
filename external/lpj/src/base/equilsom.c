/***************************************************************************/
/**                                                                       **/
/**                     e  q  u  i  l  s  o  m  . c                       **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 20.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

/*
 *  DESCRIPTION
 *
 *  Analytically solves differential flux equations for fast and slow SOM 
 *  pools assuming annual litter inputs close to long term equilibrium
 *
 */

void equilsom(Cell *cell) 
{
  int s;
  Stand *stand;
  
  foreachstand(stand,s,cell->standlist)
    equilsoil(&stand->soil);
} /* of 'equilsom' */
