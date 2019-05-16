/***************************************************************************/
/**                                                                       **/
/**                        t  r  e  e  .  h                               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 21.03.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef TREE_H /* Already included? */
#define TREE_H

/* Definition of constants */

#define allom1 100.0
#define allom2 40.0
#define allom3 0.67	/*0.5*/
#define allom4 0.3
#define reinickerp 1.6
#define k_latosa 4e3
#define wooddens 2e5
#define FPC_TREE_MAX 0.95

/* Declaration of datatypes */

typedef struct
{
  Real leaf,sapwood,root;
}Treephys;

typedef struct
{
  Real leaf,sapwood,heartwood,root,debt;
}Treephys2;

typedef struct
{
  int leaftype;			/*par16*/
  Treephys turnover;		/*turnover period (years) (9,11,12)*/
  Treephys cn_ratio;		/*C:N mass ratio (13-15) */
  Real crownarea_max;		/*tree maximum crown area (m2) (20)*/
  Treephys2 sapl;               /* sapling */
} Pfttreepar;

typedef struct
{
  /*int leafondays,leafoffdays;
  Bool leafon;*/
  Real height;
  Real crownarea;
  Real gddtw;
  Real aphen_raingreen;
  Treephys2 ind;
} Pfttree;

/* Declaration of functions */

extern void new_tree(Pft *);
extern Pft *newpftage(Pftpar *,int);
extern Real npp_tree(Pft *,Real,Real,Real,Real);
extern void mortality_tree(Litter *,Pft *,Real,Real,Real *);
extern Real fpc_tree(Pft *);
extern void litter_update_tree(Litter *,const Pft *,Real);
extern void allometry_tree(Pft *pft);
extern void allocation_tree(Litter *,Pft *,Real *);
extern Real lai_tree(const Pft *);
extern Real turnover_tree(Litter *,Pft *);
extern Real phenology_tree(Pft *,Real,Real,int);
extern Bool fwrite_tree(FILE *,const Pft *);
extern void fprint_tree(FILE *,const Pft *);
extern void fread_tree(FILE *,Pft *,Bool);
extern Bool fscanpft_tree(FILE *,Pftpar *,const char *);
extern Bool isneg_tree(const Pft *,Real);
extern Real establishment_tree(Pft *,Real,Real,int);
extern void reproduction_tree(Litter *,Pft *);
extern void init_tree(Pft *);
extern Real fire_tree(Pft *,Real);
extern Real vegc_sum_tree(const Pft *);
extern void free_tree(Pft *);
extern void light_tree(Litter *,Pft *,Real);
extern void adjust_tree(Litter *,Pft *,Real);

/* Definitions of macros */

#define istree(pft) (getpftpar(pft,type)==TREE)
#define israingreen(pft) getpftpar(pft,phenology)==RAINGREEN
#define fprinttreephys2(file,phys,nind) fprintf(file,"%6.2f %6.2f %6.2f %6.2f %6.2f",phys.leaf*nind,phys.sapwood*nind,phys.heartwood*nind,phys.root*nind,phys.debt*nind)
#define phys_sum_tree(ind) (ind.leaf+ind.root+ind.heartwood+ind.sapwood)

#endif
