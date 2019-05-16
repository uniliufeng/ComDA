#!/bin/sh
# IMPLEMENT (cd $pft && $(MAKE)) and (cd $pft && $(MAKE) clean) in src/Makefile !!!!!

########################################################################
##                                                                    ##
##               n  e  w  p  f  t  .  s  h                            ##
##                                                                    ##
##   Shell script to create new PFT templates                         ##
##                                                                    ##
##   written by Werner von Bloh                                       ##
##   Potsdam Institute for Climate Impact Research                    ##
##   P.O. Box 60 12 03                                                ##
##   14412 Potsdam/Germany                                            ##
##                                                                    ##
##   Last change: 29.11.2004                                          ##
##                                                                    ##
########################################################################

function header {

cat >src/$1/$2 << EOF
/***************************************************************************/
/**                                                                       **/
EOF
printf "/**            %16s                                           **/\n" $2 >>src/$1/$2
cat >>src/$1/$2 << EOF
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
EOF
printf "/**     Last change: %s                                           **/\n" $(date +"%d.%m.%Y") >> src/$1/$2
cat >>src/$1/$2 << EOF
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"
#include "$1.h"

EOF
}

pft=$1
mkdir -p src/$pft
cat >src/$pft/Makefile << EOF
########################################################################
##                                                                    ##
##               M  a  k  e  f  i  l  e                               ##
##                                                                    ##
EOF
printf "##   Makefile for C implementation of %-31s ##\n" "$pft part of LPJ" >>src/$pft/Makefile
cat >>src/$pft/Makefile << EOF
##                                                                    ##
##   written by Werner von Bloh                                       ##
##   Potsdam Institute for Climate Impact Research                    ##
##   P.O. Box 60 12 03                                                ##
##   14412 Potsdam/Germany                                            ##
##                                                                    ##
EOF
printf "##   Last change: %s                                          ##\n" $(date +"%d.%m.%Y") >> src/$pft/Makefile
cat >>src/$pft/Makefile << EOF
##                                                                    ##
########################################################################

include ../../Makefile.inc

OBJ	= lai_$pft.\$O turnover_$pft.\$O npp_$pft.\$O phenology_$pft.\$O\\
          reproduction_$pft.\$O mortality_$pft.\$O init_$pft.\$O\\
          allocation_$pft.\$O isneg_$pft.\$O light_$pft.\$O\\
          litter_update_$pft.\$O fpc_$pft.\$O allometry_$pft.\$O\\
          new_$pft.\$O fwrite_$pft.\$O fscanpft_$pft.\$O  fprint_$pft.\$O \\
          fread_$pft.\$O establishment_$pft.\$O fire_$pft.\$O free_$pft.\$O\\
          vegc_sum_$pft.\$O adjust_$pft.\$O

INC     = ../../include
LIB	= ../../lib/lib$pft.\$A

HDRS    = \$(INC)/buffer.h \$(INC)/coord.h \$(INC)/lpj.h \$(INC)/pftlist.h\\
          \$(INC)/soil.h \$(INC)/climate.h \$(INC)/date.h \$(INC)/pft.h\\
          \$(INC)/pftpar.h \$(INC)/types.h \$(INC)/$pft.h\\
          \$(INC)/errmsg.h \$(INC)/numeric.h\\
          \$(INC)/conf.h \$(INC)/swap.h \$(INC)/soilpar.h \$(INC)/stand.h\\
          \$(INC)/list.h \$(INC)/cell.h  \$(INC)/units.h \\
          \$(INC)/config.h

\$(LIB): \$(OBJ)
	\$(AR) \$(ARFLAGS)\$(LIB) \$(OBJ)

\$(OBJ): \$(HDRS)

.c.\$O: 
	\$(CC) \$(CFLAGS) -I\$(INC) -c \$*.c
clean:
	\$(RM) \$(RMFLAGS) \$(OBJ) \$(LIB)
EOF

cat >include/$pft.h << EOF
/***************************************************************************/
/**                                                                       **/
EOF
printf "/**            %16s                                           **/\n" $pft.h >> include/$pft.h 
cat >>include/$pft.h << EOF
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
EOF
printf "/**     Last change: %s                                           **/\n" $(date +"%d.%m.%Y") >> include/$pft.h
cat >>include/$pft.h << EOF
/**                                                                       **/
/***************************************************************************/

#ifndef ${pft}_H /* Already included? */
#define ${pft}_H

/* Definition of constants */

/* Declaration of datatypes */

typedef struct
{
} Pft$pft;

/* Declaration of functions */

extern void new_$pft(Pft *);
extern Real npp_$pft(Pft *,Real,Real,Real,Real);
extern void mortality_$pft(Litter *,Pft *,Real,Real);
extern Real fpc_$pft(Pft *);
extern void litter_update_$pft(Litter *,const Pft *,Real);
extern void allometry_$pft(Pft *pft);
extern void allocation_$pft(Litter *,Pft *,Real *);
extern Real lai_$pft(const Pft *);
extern Real turnover_$pft(Litter *,Pft *);
extern Real phenology_$pft(Pft *,Real,int,Real,int);
extern Bool fwrite_$pft(FILE *,const Pft *);
extern void fprint_$pft(FILE *,const Pft *);
extern void fread_$pft(FILE *,Pft *,Bool);
extern Bool fscanpft_$pft(FILE *,Pftpar *,const char *);
extern Bool isneg_$pft(const Pft *);
extern Real establishment_$pft(Pft *,Real,Real,int);
extern void reproduction_$pft(Litter *,Pft *);
extern void init_$pft(Pft *);
extern Real fire_$pft(Pft *,Real);
extern Real vegc_sum_$pft(const Pft *);
extern void free_$pft(Pft *);
extern void light_$pft(Litter *,Pft *,Real);
extern Real adjust_$pft(Litter *,Pft *,Real);

/* Definitions of macros */


#endif
EOF
header $pft new_$pft.c
cat >>src/$pft/new_$pft.c  << EOF
void new_$pft(Pft *pft)
{
  Pft$pft *$pft;
  $pft=new(Pft$pft);
  pft->data=$pft; 
} /* of 'new_$pft' */
EOF
header $pft npp_$pft.c
cat >>src/$pft/npp_$pft.c  << EOF
Real npp_$pft(Pft *pft,Real phen,Real gtemp_air,Real gtemp_soil,Real assim)
{
} /* of 'npp_$pft' */
EOF
header $pft mortality_$pft.c
cat >>src/$pft/mortality_$pft.c  << EOF
void mortality_$pft(Litter *litter,Pft *pft,Real turnover_ind,Real mtemp_max)
{
} /* of 'mortality_$pft' */
EOF
header $pft fpc_$pft.c
cat >>src/$pft/fpc_$pft.c  << EOF
Real fpc_$pft(Pft *pft)
{
} /* of 'fpc_$pft' */
EOF
header $pft litter_update_$pft.c
cat >>src/$pft/litter_update_$pft.c  << EOF
void litter_update_$pft(Litter *litter,const Pft *pft,Real frac)
{
} /* of 'litter_update_$pft' */
EOF
header $pft allometry_$pft.c
cat >>src/$pft/allometry_$pft.c  << EOF
void allometry_$pft(Pft *pft)
{
} /* of 'allometry_$pft' */
EOF
header $pft allocation_$pft.c
cat >>src/$pft/allocation_$pft.c  << EOF
void allocation_$pft(Litter *litter,Pft *pft,Real *fpc_inc)
{
} /* of 'allocation_$pft' */
EOF
header $pft lai_$pft.c
cat >>src/$pft/lai_$pft.c  << EOF
Real lai_$pft(const Pft *pft)
{
} /* of 'lai_$pft' */
EOF
header $pft turnover_$pft.c
cat >>src/$pft/turnover_$pft.c  << EOF
Real turnover_$pft(Litter *litter,const Pft *pft)
{
} /* of 'turnover_$pft' */
EOF
header $pft phenology_$pft.c
cat >>src/$pft/phenology_$pft.c  << EOF
Real phenology_$pft(Pft *pft,Real temp,int gdd5,Real lat,int day)
{
} /* of 'phenology_$pft' */
EOF
header $pft fwrite_$pft.c
cat >>src/$pft/fwrite_$pft.c  << EOF
Bool fwrite_$pft(FILE *file,const Pft *pft)
{
  return (fwrite(pft->data,sizeof(Pft$pft),1,file)!=1);
} /* of 'fwrite_$pft' */
EOF
header $pft fprint_$pft.c
cat >>src/$pft/fprint_$pft.c  << EOF
void fprint_$pft(FILE *file,const Pft *pft)
{
} /* of 'fprint_$pft' */
EOF
header $pft fread_$pft.c
cat >>src/$pft/fread_$pft.c  << EOF
void fread_$pft(FILE *file,Pft *pft,Bool swap)
{
} /* of 'fread_$pft' */
EOF
header $pft fscanpft_$pft.c
cat >>src/$pft/fscanpft_$pft.c  << EOF
Bool fscanpft_$pft(FILE *file,Pftpar *pft,const char *filename)
{
  pft->newpft=new_$pft;
  pft->turnover=turnover_$pft;
  pft->npp=npp_$pft;
  pft->leaf_phenology=phenology_$pft;
  pft->fwrite=fwrite_$pft;
  pft->fprint=fprint_$pft;
  pft->fread=fread_$pft;
  pft->litter_update=litter_update_$pft;
  pft->allocation=allocation_$pft;
  pft->isneg=isneg_$pft;
  pft->establishment=establishment_$pft;
  pft->reproduction=reproduction_$pft;
  pft->mortality=mortality_$pft;
  pft->init=init_$pft;
  pft->fire=fire_$pft;
  pft->lai=lai_$pft;
  pft->free=free_$pft;
  pft->vegc_sum=vegc_sum_$pft;
  pft->light=light_$pft;
  pft->adjust=adjust_$pft;
  return FALSE;
} /* of 'fscanpft_$pft' */
EOF
header $pft isneg_$pft.c
cat >>src/$pft/isneg_$pft.c  << EOF
Bool isneg_$pft(const Pft *pft)
{
} /* of 'isneg_$pft' */
EOF
header $pft establishment_$pft.c
cat >>src/$pft/establishment_$pft.c  << EOF
Real establishment_$pft(Pft *pft,Real fpc_total,Real fpc_tree,int n_est)
{
} /* of 'establishment_$pft' */
EOF
header $pft reproduction_$pft.c
cat >>src/$pft/reproduction_$pft.c  << EOF
void reproduction_$pft(Litter *litter,Pft *pft)
{
} /* of 'reproduction_$pft' */
EOF
header $pft init_$pft.c
cat >>src/$pft/init_$pft.c  << EOF
void init_$pft(Pft *pft)
{
} /* of 'init_$pft' */
EOF
header $pft fire_$pft.c
cat >>src/$pft/fire_$pft.c  << EOF
Real fire_$pft(Pft *pft,Real fireprob)
{
} /* of 'fire_$pft' */
EOF
header $pft vegc_sum_$pft.c
cat >>src/$pft/vegc_sum_$pft.c  << EOF
Real vegc_sum_$pft(const Pft *pft)
{
} /* of 'vegc_sum_$pft' */
EOF
header $pft free_$pft.c
cat >>src/$pft/free_$pft.c  << EOF
void free_$pft(Pft *pft)
{
  free(pft->data);
} /* of 'free_$pft' */
EOF
header $pft light_$pft.c
cat >>src/$pft/light_$pft.c  << EOF
void light_$pft(Litter *litter,Pft *pft,Real excess)
{
} /* of 'light_$pft' */
EOF
header $pft adjust_$pft.c
cat >>src/$pft/adjust_$pft.c  << EOF
Real adjust_$pft(Litter *litter,Pft *pft,Real fpc_$pft)
{
} /* of 'adjust_$pft' */
EOF
