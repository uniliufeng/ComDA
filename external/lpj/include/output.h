/***************************************************************************/
/**                                                                       **/
/**                        o  u  t  p  u  t  .  h                         **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 01.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef OUTPUT_H
#define OUTPUT_H

/* Definition of datatypes */

typedef struct
{
  Real mnpp;
  Real mrh;
  Real mtransp;
  Real mrunoff;
  Real mevap;
  Real minterc;
  Real mswc[NSOILLAYER];
  Real firec;
  Real firef;
  Real flux_estab;
} Output; 

/* Declaration of functions */

extern FILE **fopenoutput(Config);
extern void fcloseoutput(FILE **,int);
extern void initoutput(Output *);

#endif
