/***************************************************************************/
/**                                                                       **/
/**                     s  w  a  p  .  h                                  **/
/**                                                                       **/
/**     Converts big endian into little endian and vice versa             **/
/**     Needed for reading data under IA32                                **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  14.10.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#ifndef SWAP_H /* Already included? */
#define SWAP_H

/* Declaration of functions */

extern short int swapshort(short int);
extern int swapint(int);
extern int freadint(int *,int,Bool,FILE *);
extern int freadreal(Real *,int,Bool,FILE *);

/* Definitions of macros */

#define freadint1(data,swap,file) freadint(data,1,swap,file)
#define freadreal1(data,swap,file) freadreal(data,1,swap,file)

#endif
