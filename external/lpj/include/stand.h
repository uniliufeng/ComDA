/***************************************************************************/
/**                                                                       **/
/**                    s  t  a  n  d  .  h                                **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 27.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef STAND_H  /* Already included? */
#define STAND_H

/* Definition of datatypes */

typedef struct
{
  Pftlist  pftlist;
  Soil soil;
  Real fire_sum;
  Real totc2;
  Real frac;
}Stand;

typedef List *Standlist;

/* Declaration of functions */

extern Bool fwritestand(FILE *,const Stand *,Bool);
extern void fprintstand(FILE *,const Stand *);
extern int fwritestandlist(FILE *,const Standlist,Bool);
extern void fprintstandlist(FILE *,const Standlist);
extern Stand *freadstand(FILE *,const Pftpar[],const Soilpar *,Bool);
extern Standlist freadstandlist(FILE *,const Pftpar [],const Soilpar *,Bool);

/* Definition of macros */

#define getstand(list,index) ((Stand *)getlistitem(list,index))
#define foreachstand(stand,i,list) for(i=0;i<getlistlen(list) && (stand=getstand(list,i));i++)
#define printstandlist(standlist) fprintstandlist(stdout,standlist)

#endif
