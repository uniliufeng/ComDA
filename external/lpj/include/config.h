/***************************************************************************/
/**                                                                       **/
/**                   c  o  n  f  i  g  .  h                              **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 14.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef CONFIG_H /* Already included? */
#define CONFIG_H

/* Definition of datatypes */

typedef struct
{
  char *pftpar_filename;
  char *soilpar_filename;
  char *coord_filename;
  char *soil_filename;
  char *temp_filename;
  char *prec_filename;
  char *cloud_filename;
  char *wet_filename;
  char *co2_filename;
  char *restart_filename;
  char *write_restart_filename;
  char **out_filename;
  int ngridcell,startgrid,nspinup,lastyear,totalgridcell,firstyear;
  int n_out;
  int seed;
  Coord resolution;
} Config;

/* Declaration of functions */

extern Bool fscanconfig(Config *,int *,char***);

#endif
