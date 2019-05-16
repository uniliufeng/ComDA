/***************************************************************************/
/**                                                                       **/
/**                      c  e  l  l  .  h                                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     LPJ header file contains all necessary includes                   **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef CELL_H /* Already included? */
#define CELL_H

/* Definition of datatypes */

typedef struct
{
  Coord coord;
  Standlist standlist;
  Climbuf climbuf;
  Real *gdd;
  Real aprec;
  Bool skip;
  Output output;
} Cell;

/* Declaration of functions */

extern void freegrid(Cell *,int);
extern void update_daily(Cell *,Real,Real,Real,Real,int);
extern Real update_annual(FILE **,Cell *,const Pftpar [],int,int,int,Config);
extern Real update_monthly(FILE **,Cell *,Real,int,int,int,Config);
extern void init_annual(Cell *,int);
extern int fwritecell(FILE *,int [],const Cell [],int,int,Bool);
extern void fprintcell(FILE *,const Cell [],int);
extern Bool freadcell(FILE *,Cell *,const Pftpar[],int,const Soilpar *,Bool);
extern int writecoords(FILE *,const Cell [],int);
extern int iterate(FILE **,Cell *,Climate *,const Pftpar *,int,int,Config);
extern void iterateyear(FILE **,Cell *,Climate *,Real,const Pftpar *,int,
                        int,Config,int);
extern int fwriteoutput_annual(FILE **,const Cell *,int,int);
extern int fwriteoutput_monthly(FILE **,const Cell *,int,int);
extern void equilsom(Cell *);

#define printcell(cell,n) fprintcell(stdout,cell,n)

#endif
