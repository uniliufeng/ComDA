/***************************************************************************/
/**                                                                       **/
/**                     h  e  a  d  e  r  .  h                            **/
/**                                                                       **/
/**     Reading/Writing file header for  LPJ related files. Detects       **/
/**     whether byte order has to be changed                              **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  29.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#ifndef HEADER_H /* Already included? */
#define HEADER_H

/*Definition of constants */

#define RESTART_HEADER "LPJRESTART"
#define RESTART_VERSION 4

/* Definition of datatypes */

typedef struct
{
  int year,firstcell,ncell;
} Restartheader;

/* Declaration of functions */

extern Bool fwriteheader(FILE *,const char *,int);
extern Bool freadheader(FILE *,int *,const char *,int);

#endif
