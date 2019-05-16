/***************************************************************************/
/**                                                                       **/
/**                     p  f  t  .  h                                     **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 31.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef PFT_H /* Already included? */
#define PFT_H

#include "pftpar.h"


/* Definitions of datatypes */

typedef struct
{
  Real low,high;
} Limit;


typedef struct Pft
{
  const struct  Pftpar
  {
    int id,type;	      	/* type --> whether CROP or TREE or GRASS*/
    char *name;              	/* Pft name */
    Real rootdist[NSOILLAYER];	/* fraction of roots in upper soil layer par1*/
    Real minwscal;   	    	/* water scalar value at which leaves shed by 
                              	   drought deciduous PFT par3 */
    Real gmin;	            	/* canopy conductance component (4) */
    Real respcoeff;           	/* maintenance respiration coefficient (5) */
    Real nmax;			/* maximum foliar N content (mg/g) (7) */
    Real resist;		/* fire resistance index (8) */
    Real longivity;		/* leaf longivity (10) */
    Real lmro_ratio;            /* leaf to root ratio under non-water stressed 
                                   conditions (18) */
    Real ramp;              	/* number of GDDs to attain full leaf cover 
                                   (par19) */
    Real gdd5min;         	/* PFT-specific minimum GDD(30) */
    Real twmax;                 /* (31) */
    Real gddbase;		/* GDD base  (33) */
    Real min_temprange;         /*  (34) */
    Real emax;			/* emax (mm/day)(35) */
    Real intc;			/* Interception storage parameter (36) */
    Real lai_sapl;
    int phenology;              /* par17 */
    int path;			/* par2 */
    Real sla;
    Limit temp_co2;		/* temperature limit for CO2 uptake (24,27) */
    Limit temp_photos;    	/* range of temperature optimum for 
                                   photosynthesis(25,26) */
    Limit temp;        		/* bioclimatic limits (28,29) */
    void *data;                 /* pointer for PFT specific extensions */

    /* list of pointers for PFT specific functions */

    void (*newpft)(struct Pft *); 
    void (*init)(struct Pft *);   
    Real (*turnover)(Litter *,struct Pft *);
    Real (*npp)(struct Pft*,Real,Real,Real,Real);
    Real (*fpc)(struct Pft*);
    Real (*leaf_phenology)(struct Pft *,Real,Real,int);
    Bool (*fwrite)(FILE *,const struct Pft *);
    void (*fread)(FILE *,struct Pft *,Bool);
    void (*fprint)(FILE *,const struct Pft *);
    void (*litter_update)(Litter *,const struct Pft *,Real);
    void (*allocation)(Litter *,struct Pft *,Real *);
    Bool (*isneg)(const struct Pft *, Real);
    Real (*establishment)(struct Pft *,Real,Real,int);
    void (*reproduction)(Litter *,struct Pft *);
    void (*mortality)(Litter *,struct Pft *,Real,Real,Real *);
    Real (*fire)(struct Pft *,Real);
    Real (*lai)(const struct Pft *);
    void (*adjust)(Litter *,struct Pft *,Real);
    void (*free)(struct Pft *);
    void (*light)(Litter *,struct Pft *,Real);
    Real (*vegc_sum)(const struct Pft *);
  } *par;                /* Pft parameters */
  Real fpc;              /* foliar projective cover (FPC) under full leaf 
                            cover as fraction of modelled area */
  Real nind;
  Real gdd;
  Real bm_inc;
  Real wscal;
  Real wscal_mean;
  Real phen,aphen;
  Buffer *gddbuf;
  void *data;            /* pointer for PFT specific extensions */
} Pft;

typedef struct Pftpar Pftpar;


typedef Bool (*Fscanpftparfcn)(FILE *,Pftpar *,const char *);

/* Declaration of functions */

extern Pft *newpft(const Pftpar *);
extern void initpft(Pft *);
extern void freepft(Pft *);
extern int fscanpftpar(Pftpar **,const char *,const Fscanpftparfcn [],int);
extern Bool establish(Real,const Pftpar *,const Climbuf *);
extern Bool survive(const Pftpar *,const Climbuf *);
extern Real temp_stress(const Pftpar *,Real,Real);
extern Real photosynthesis(Real *,Real *,int,Real,Real,Real,Real ,Real,Real,Real);
extern Bool killpft(Litter *,Pft *,const Climbuf *);
extern Real interception(Real *,const Pft *,Real,Real);
extern void initgdd(Real [],int);
extern void updategdd(Real [],const Pftpar [],int,Real);
extern Real gp(Pft *,Real,Real,Real,Real,Real *,Real *);
extern Real water_stressed(Pft *,Real [NSOILLAYER],const Real [NSOILLAYER],
                           Real,Real,Real,Real *,Real *,Real,Real,Real,Real,Real);
extern Bool fwritepft(FILE *,const Pft *,Bool);
extern void fprintpft(FILE *,const Pft *);
extern Pft *freadpft(FILE *,const Pftpar[],Bool);
extern void nomortality(Litter *,Pft *,Real,Real,Real *);
extern void noinit(Pft *);
extern Real nofire(Pft *,Real);
extern void noadjust(Litter *,Pft *,Real);

/* Definition of macros */

#define isphoto(tstress) (tstress>=1e-2)
#define getpftpar(pft,val) (pft)->par->val
#define newgdd(npft) newvec(Real,npft)
#define fpc(pft)
#define npp(pft,phen,gtemp_air,gtemp_soil,assim) pft->par->npp(pft,phen,gtemp_air,gtemp_soil,assim)
#define leaf_phenology(pft,temp,lat,day) pft->par->leaf_phenology(pft,temp,lat,day)
#define allocation(litter,pft,fpc_inc) pft->par->allocation(litter,pft,fpc_inc)
#define reproduction(litter,pft) pft->par->reproduction(litter,pft)
#define mortality(litter,pft,turnover_ind,mtemp_max,fpc_inc) pft->par->mortality(litter,pft,turnover_ind,mtemp_max,fpc_inc)
#define isneg(pft,bm_inc) pft->par->isneg(pft,bm_inc)
#define litter_update(litter,pft,frac) pft->par->litter_update(litter,pft,frac)
#define turnover(litter,pft) pft->par->turnover(litter,pft)
#define fire(pft,fireprob) pft->par->fire(pft,fireprob)
#define lai(pft) pft->par->lai(pft)
#define vegc_sum(pft) pft->par->vegc_sum(pft)
#define adjust(litter,pft,fpc) pft->par->adjust(litter,pft,fpc)
#define establishment(pft,fpc_total,fpc,n_est) pft->par->establishment(pft,fpc_total,fpc,n_est);

#endif
