/***************************************************************************/
/**                                                                       **/
/**                     t  y  p  e  s  .  h                               **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran version         **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 04.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef TYPES_H /* Already included? */
#define TYPES_H

/* Definition of constants */

#ifndef TRUE    /* Check whether TRUE or FALSE are already defined */
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define STRING_LEN 255 /* Maximum length of default strings  */

/* Definition of datatypes */

typedef double Real; /* Default floating point type in LPJ */

typedef int Bool;

typedef char String[STRING_LEN+1];

/* Declaration of functions */

extern void fail(const char *,...);
extern Bool fscanreal(FILE *,Real *);
extern Bool fscanstring(FILE *,char *);

/* Definitions of macros */

#define newvec(type,size) (type *)malloc(sizeof(type)*(size))
#define new(type) (type *)malloc(sizeof(type))
#ifndef min
#define min(a,b) ((a<b) ?  a : b)
#endif
#ifndef max
#define max(a,b) ((a>b) ?  a : b)
#endif
#define fwrite1(var,size,file) if(fwrite(var,size,1,file)!=1) return TRUE
#define fwriten(var,size,n,file) if(fwrite(var,size,n,file)!=n) return TRUE
#define setminthreshold(a,threshold) if(a<threshold) a=0

#endif /* of TYPES_H */
