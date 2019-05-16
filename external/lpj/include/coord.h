/***************************************************************************/
/**                                                                       **/
/**                    c  o  o  r  d  .  h                                **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 05.10.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef COORD_H /* Already included? */
#define COORD_H

/* Definition of datatypes */

typedef struct
{
  Real lon,lat;  /* longitude, latitude in degrees */
  Real area;     /* cell area (m^2) */
} Coord;
typedef struct
{
  FILE *file;
  int n;
  Bool swap;
} Coordfile;

/* Declaration of functions */

extern Coordfile *opencoord(const char *);
extern int seekcoord(Coordfile *,int);
extern Bool readcoord(Coordfile *,Coord *,Coord);
extern void closecoord(Coordfile *);
extern Bool writecoord(FILE *,Coord);
extern Real cellarea(Coord,Coord);

/* Definition of macros */

#define fprintcoord(file,coord) fprintf(file,"%5.2f %5.2f",coord.lon,coord.lat)
#define fscancoord(file,coord) fscanf(file,"%lg %lg",coord.lon,coord.lat)
#define numcoord(coordfile) coordfile->n

#endif
