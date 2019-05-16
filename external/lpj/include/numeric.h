/***************************************************************************/
/**                                                                       **/
/**                   n  u  m  e  r  i  c  .  h                           **/
/**                                                                       **/
/**     Header for numerical utility routines                             **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 18.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef NUMERIC_H /* Already included? */
#define NUMERIC_H

/* Definition of datatypes */

typedef Real (*Bisectfcn)(Real,void *);

/* Declaration of functions */

extern Real bisect(Real (*)(Real,void *),Real,Real,void *,Real,Real,int); /* find zero */
extern Real leftmostzero(Real (*)(Real,void *),Real,Real,void *,Real,Real,int); /* find leftmost zero */
extern void linreg(Real *,Real *,const Real[],int); /* linear regression */
extern void setseed(int); /* set seed of random number generator */
extern int getseed(void); /* get seed of random number generator */
extern Real randfrac(void); /* random number generator */
extern void petpar(Real *,Real *,Real *,Real,Real,Real,Real);

#endif
