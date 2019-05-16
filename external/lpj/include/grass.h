/***************************************************************************/
/**                                                                       **/
/**                   g  r  a  s  s  .  h                                 **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef GRASS_H /* Already included? */
#define GRASS_H


/* Definition of datatypes */

typedef struct
{
  Real leaf,root;
}Grassphys;

typedef struct
{
  Grassphys turnover;		/*turnover period (years) (9,11,12)*/
  Grassphys cn_ratio;		/*C:N mass ratio (13-15) */
  Grassphys sapl;               /* sapling */
} Pftgrasspar;

typedef struct
{
  Grassphys ind;
} Pftgrass;

/* Declaration of functions */

extern void new_grass(Pft *);
extern Real npp_grass(Pft *,Real,Real,Real,Real);
extern Real fpc_grass(Pft *);
extern void litter_update_grass(Litter *,const Pft*,Real);
extern void allocation_grass(Litter *,Pft *,Real *);
extern Real lai_grass(const Pft *);
extern Real turnover_grass(Litter *,Pft *);
extern Real phenology_grass(Pft *,Real,Real,int);
extern Bool fwrite_grass(FILE *,const Pft *);
extern void fprint_grass(FILE *,const Pft *);
extern void fread_grass(FILE *,Pft *,Bool);
extern Bool fscanpft_grass(FILE *,Pftpar *,const char *);
extern Bool isneg_grass(const Pft *, Real);
extern Real establishment_grass(Pft *,Real,Real,int);
extern Real vegc_sum_grass(const Pft *);
extern void reproduction_grass(Litter *,Pft *);
extern void init_grass(Pft *);
extern void free_grass(Pft *);
extern void light_grass(Litter *,Pft *,Real);

/* Definition of macros */

#define isgrass(pft) (getpftpar(pft,type)==GRASS)
#define fprintgrassphys(file,phys) fprintf(file,"%6.2f %6.2f",phys.leaf,phys.root)
#define phys_sum_grass(ind) (ind.leaf+ind.root)

#endif
