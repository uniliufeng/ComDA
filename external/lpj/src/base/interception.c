/***************************************************************************/
/**                                                                       **/
/**                  i  n  t  e  r  c  e  p  t  i  o  n  .  c             **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 30.05.2005                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

#define PT 1.32	 /* Priestley-Taylor coefficient */

Real interception(Real *wet,const Pft *pft,Real pet,Real rain)   
{
  Real int_store;
  if(pet<0.0001)
  {
    *wet=0;
    return 0; 
  }
  int_store=pft->par->intc*lai(pft)*pft->phen;
  if(int_store>0.9999)
    int_store=0.9999;
  *wet=int_store*rain/(pet*PT);
  if(*wet>0.9999)
    *wet=0.9999;
  return pet*PT*(*wet)*pft->fpc;
} /* of interception */
