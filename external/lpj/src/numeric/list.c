/***************************************************************************/
/**                                                                       **/
/**                       l  i  s  t  .  c                                **/
/**                                                                       **/
/**     C implementation of a resizeable array                            **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change:  10.09.2004                                          **/
/**                                                                       **/
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "types.h"
#include "list.h"

List *newlist(void)
{
  List *list;
  list=(List *)malloc(sizeof(List));
  list->n=0;
  list->data=NULL;
  return list;
} /* of 'newlist' */

int addlistitem(List *list,void *item)
{
  list->data=(void **)realloc(list->data,sizeof(void *)*(list->n+1));
  list->data[list->n++]=item;
  return list->n;
} /* of 'addlistitem' */

int dellistitem(List *list,int i)
{
  /* does not check for empty list if SAFE not defined */
#ifdef SAFE
  if(isempty(list))
    fail("list is empty!\n");
#endif
  list->n--;
  list->data[i]=list->data[list->n];
  if(isempty(list))
  {
    free(list->data);
    list->data=NULL;
  }
  else
    list->data=(void **)realloc(list->data,sizeof(void *)*list->n);
  return list->n;
} /* of 'dellistitem' */

void freelist(List *list)
{
  if(!isempty(list))
    free(list->data);
  free(list);
} /* of 'freelist' */
