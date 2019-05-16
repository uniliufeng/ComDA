/***************************************************************************/
/**                                                                       **/
/**                         c  r  o  p  .  h                              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.05.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef CROP_H /* Already included? */
#define CROP_H

/* Definition of datatypes */

typedef struct
{
  Real frac;
} Pftcrop;

typedef struct
{
  Real frac;
}  Pftcroppar;

/* Declaration of functions */

void npp_crop(Pft *,Real,Real);

/* Definition of macros */

#define iscrop(pft) (getpftpar(pft,type)==CROP)

#endif
