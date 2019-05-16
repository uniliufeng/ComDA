/***************************************************************************/
/**                                                                       **/
/**                  p  f  t  l  i  s  t  .  c                            **/
/**                                                                       **/
/**     C implementation of LPJ, derived from the Fortran/C++ version     **/
/**                                                                       **/
/**     written by Werner von Bloh, Sibyll Schaphoff                      **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include "lpj.h"

int delpft(Pftlist pftlist, /* Pft list */
           int index        /* index of entry to delete */
          )                 /* returns updated number of entries in pft list */
{
  freepft(getlistitem(pftlist,index));
  return dellistitem(pftlist,index);
} /* of 'delpft ' */

int fwritepftlist(FILE *file,const Pftlist pftlist,Bool full)
{
  const Pft *pft;
  int p;
  fwrite(&pftlist->n,sizeof(int),1,file);
  foreachpft(pft,p,pftlist)
    if(fwritepft(file,pft,full))
      break;
  return p;
} /* of 'fwritepftlist' */

void fprintpftlist(FILE *file,const Pftlist pftlist)
{
  const Pft *pft;
  int p;
  fprintf(file,"Number of Pfts: %d\n",pftlist->n);
  foreachpft(pft,p,pftlist)
    fprintpft(file,pft);
} /* of 'fprintpftlist' */

Pftlist freadpftlist(FILE *file,const Pftpar pftpar[],Bool swap)
{
  int p,n;
  Pftlist pftlist;
  if(freadint(&n,1,swap,file)!=1)
    return NULL;
  pftlist=newpftlist();
  for(p=0;p<n;p++)
    addlistitem(pftlist,freadpft(file,pftpar,swap));
  return pftlist;
} /* of 'freadpftlist' */

void freepftlist(Pftlist pftlist /* Pft list */
                )
{
  if(!isempty(pftlist))
    while(delpft(pftlist,0)>0); 
  freelist(pftlist);
} /* of 'freepftlist' */
