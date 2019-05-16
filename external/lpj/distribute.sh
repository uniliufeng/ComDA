#!/bin/sh

#############################################################################
##                                                                         ##
##                d  i  s  t  r  i  b  u  t  e  .  s  h                    ##
##                                                                         ##
## sh script to create LPJ config files for distribiuted grid              ##
##                                                                         ##
## Usage:  distribute.sh jobnumber [gridnumber]                            ##
##                                                                         ##
## Last change: 14.10.2004                                                 ##
##                                                                         ##
#############################################################################

if [ $# -lt 1 ] 
then
  echo >&2 usage: distribute.sh jobnumber [gridnumber]
  exit 1
fi
if [ $# -gt 1 ]
then
  n=$2
else
  n=67420
fi
p=$1
if [ $p -lt 1 ]
then
  echo >&2 Error: number of jobs less than one.
  exit 1
fi
index=0
outdir=output
while (($index < $p))
  do
  start=$((  n*index / p ))
  stop=$(( n*(index+1)/p-1 ))
  cat <<EOF >lpj$index.cmd
#!/bin/ksh
# @ executable = lpj_test 
# @ arguments = -DISRANDOM lpj$index.conf
# @ input = /dev/null
# @ output = lpj$index.out
# @ error = lpj$index.err
# @ resources = ConsumableCpus(1)
# @ initialdir =/scratch/01/sibyll/new_code/LPJ-C-0.09.28
# @ notify_user = sibylls@pik-potsdam.de
# @ class = short
# @ group = bis
# @ notification = complete
# @ checkpoint = no
# @ restart = no
# @ requirements = (Arch == "R6000") && (OpSys == "AIX51")  
# @ queue
EOF
  cat  >lpj$index.conf <<EOF
/*********************************************************************/
/*                                                                   */
/*                   l  p  j  $index  .  c  o  n  f                  */
/*                                                                   */
/* Configuration file for LPJ C Version, created by distribute.sh    */
/*                                                                   */
/* Last change: 29.09.2004                                           */
/*                                                                   */
/*********************************************************************/

#include "include/conf.h" /* include constant definitions */


/*#define ISRANDOM*/  /* random generation of daily precipitation */


par/pft.par  /* PFT Parameter file */
par/soil.par  /* Soil Parameter file */
#include "cru_1901-2003_test.conf"  /* Input files of CRU dataset */
#ifdef ISRANDOM
RANDOM_SEED
#endif
16	             /* number of outputfiles */
$outdir/grid$index.out
$outdir/fpc$index.bin
$outdir/mnpp$index.bin
$outdir/mrh$index.bin
$outdir/mtransp$index.bin
$outdir/mrunoff$index.bin
$outdir/mevap$index.bin
$outdir/minterc$index.bin
$outdir/mswc1$index.bin
$outdir/mswc2$index.bin
$outdir/firec$index.bin
$outdir/firef$index.bin
$outdir/vegc$index.bin
$outdir/soilc$index.bin
$outdir/litc$index.bin
$outdir/flux_estab$index.bin
FIRE /* fire disturbance enabled */ 
$start /* first grid cell */
$stop /* last grid cell  */

#ifndef FROM_RESTART

1000 /* spinup years */
1901 /* first year of simulation */
1901 /* last year of simulation */
NO_RESTART /* do not start from restart file */
RESTART /* create restart file: the last year of simulation=restart-year */
$outdir/restart$index.lpj /* filename of restart file */

#else

0    /* no spinup years */
1901 /* first year of simulation */
2003 /* last year of simulation */
RESTART /* start from restart file */
$outdir/restart$index.lpj
RESTART /* create restart file */
$outdir/restart_final$index.lpj   /* filename of restart file */

#endif
EOF
  llsubmit lpj$index.cmd
  let index=index+1
done
