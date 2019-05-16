/***************************************************************************/
/**                                                                       **/
/**                         s  o  i  l  .  h                              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef SOIL_H /* Already included? */
#define SOIL_H

/* Definition of constants */

#define NSOILLAYER 2 /* Number of soil layers */
#define TOPLAYER 0
#define BOTTOMLAYER (NSOILLAYER-1)
#define soil_equil_year 400
#define atmfrac 0.7
#define soilfrac (1-atmfrac)
#define fastfrac 0.98

extern  Real soildepth[NSOILLAYER];

#include "soilpar.h"

/* Definition of datatypes */

typedef struct
{
  Real fast,slow;
} Pool;

typedef struct
{
  Real ag_tree,ag_grass,bg;
} Litter;

typedef struct
{
  int type;
  char *name; /* soil name */
  Real k1,k2;
  Real whc[NSOILLAYER],whcs[NSOILLAYER];
  Real tdiff_0,tdiff_15,tdiff_100;
  Real tdiff;
} Soilpar;

typedef struct
{
  const Soilpar *par;
  Pool cpool,k_mean;
  Real w[NSOILLAYER];
  Real w_evap;
  Real alag,amp,meanw1;
  Real snowpack;
  Real decomp_litter_mean;
  Litter litter;
} Soil; 

/* Declaration of functions */

extern int fscansoilpar(Soilpar **,const char *);
extern void initsoil(Soil *soil,const Soilpar *);
extern void snow(Real *,Real *,Real *,Real);
extern Real soiltemp(const Soil *,const Climbuf *);
extern Real temp_response(Real);
extern Real littersom(Soil *,Real);
extern Real waterbalance(Soil *,Real [NSOILLAYER],Real *,Real,Real,Real);
extern void getlag(Soil *,int);
extern Bool fwritesoil(FILE *,const Soil *,Bool);
extern Bool freadsoil(FILE *,Soil *,const Soilpar *,Bool);
extern Real fire_sum(const Litter *,Real);
extern Real fire_prob(const Litter *,Real,Real *);
extern void fprintsoil(FILE *,const Soil *);
extern void equilsoil(Soil *);

/* Definition of macros */

#define getsoilpar(soil,var) (soil)->par->var
#define foreachsoillayer(l) for(l=0;l<NSOILLAYER;l++)
#define fprintpool(file,pool) fprintf(file,"%5.2f %5.2f",pool.slow,pool.fast)
#define fprintlitter(file,litter) fprintf(file,"%5.2f %5.2f %5.2f",litter.ag_tree,litter.ag_grass,litter.bg)
#define littersum(litter) (litter.bg+litter.ag_tree+litter.ag_grass)

#endif
