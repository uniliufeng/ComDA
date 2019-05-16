/***************************************************************************/
/**                                                                       **/
/**                     p  o  p  e  n  .  h                               **/
/**                                                                       **/
/**     Popen implementation for Windows/Unix                             **/
/**                                                                       **/
/**     written by Werner von Bloh                                        **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 11.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#ifndef POPEN_H
#define POPEN_H

#ifdef _WIN32
#define cpp_cmd "cl /EP "
#else
#define cpp_cmd "cpp -P "
#endif

extern FILE* pt_popen(const char *,const char *);
extern void pt_pclose(FILE *);

#endif
