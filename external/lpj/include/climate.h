/***************************************************************************/
/**                                                                       **/
/**                 c  l  i  m  a  t  e  .  h                             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  30.05.2005                                          **/
/**                                                                       **/
/***************************************************************************/

#ifndef CLIMATE_H /* Already included? */
#define CLIMATE_H

#define CLIMBUFSIZE 20 /* size of climate buffer */
#define NDAYS 31
#define CELLYEAR 1 /* order in climate data set */
#define YEARCELL 2

/* Definitions of datatypes */

typedef struct
{
  int order,firstyear,nyear,ncell;
} Climateheader;

typedef struct
{
  int n,offset,size,firstyear;
  Bool swap;
  Real *co2;
  FILE *file_temp,*file_prec,*file_cloud,*file_wet;
  Real *temp,*prec,*cloud,*wet;
} Climate;

typedef struct
{
  Real gdd5;
  Real temp[NDAYS];
  Real dval_prec[NDAYS+1];
  Real temp_min,temp_max,atemp_mean;
  Buffer *min,*max;
} Climbuf;

/* Definitions of macros */

#define getcelltemp(climate,cell) climate->temp+(cell)*NMONTH
#define getcellprec(climate,cell) climate->prec+(cell)*NMONTH
#define getcellsun(climate,cell) climate->cloud+(cell)*NMONTH
#define getcellwet(climate,cell) climate->wet+(cell)*NMONTH
#define getprec(cell,d) (cell).climbuf.dval_prec[(d)+1]
#define initgdd5(climbuf) climbuf.gdd5=0
#define updategdd5(climbuf,temp) if(temp>5) (*climbuf).gdd5++

/* Declaration of functions */
 
extern Climate *initclimate(Config);
extern Bool getclimate(Climate *,int);
extern Real getco2(const Climate *,int);
extern void freeclimate(Climate *);
extern void new_climbuf(Climbuf *);
extern void init_climbuf(Climbuf *);
extern void daily_climbuf(Climbuf *,Real);
extern void monthly_climbuf(Climbuf *,Real);
extern void annual_climbuf(Climbuf *);
extern Bool fwriteclimbuf(FILE *,const Climbuf *);
extern Bool freadclimbuf(FILE *,Climbuf *,Bool);
extern void freeclimbuf(Climbuf *);
extern void prdaily(Climbuf *,int,Real,Real);
extern Bool fwriteclimateheader(FILE *,Climateheader);
extern Bool freadclimateheader(FILE *,Climateheader *,Bool *);
extern int  climateheadersize(void);

#endif
