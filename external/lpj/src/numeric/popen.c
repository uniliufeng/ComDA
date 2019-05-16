/***************************************************************************/
/**                                                                       **/
/**                   p  o  p  e  n  .  c                                 **/
/**                                                                       **/
/**     Popen implementation for Windows 32bit                            **/
/**                                                                       **/
/**     written by Kurt Keller, Paritysoft GmbH                           **/
/**     modified by Werner von Bloh                                       **/
/**     Potsdam Institute for Climate Impact Research                     **/
/**     PO Box 60 12 03                                                   **/
/**     14412 Potsdam/Germany                                             **/
/**                                                                       **/
/**     Last change: 15.11.2004                                           **/
/**                                                                       **/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include "popen.h"

/*---------------------------------------------------------------------------
  Globals fuer die Routinen pt_popen() / pt_pclose()
  ---------------------------------------------------------------------------*/
#if defined(USE_CPP) && defined(_WIN32) 

#include <windows.h>

HANDLE my_pipein[2], my_pipeout[2], my_pipeerr[2];
char   my_popenmode = ' ';

static int my_pipe(HANDLE *readwrite)
{
  SECURITY_ATTRIBUTES sa;

  sa.nLength = sizeof(sa);          /* Laenge in Byte */
  sa.bInheritHandle = 1;            /* Descriptoren sollen vererbbar sein */
  sa.lpSecurityDescriptor = NULL;

  if (! CreatePipe (&readwrite[0],&readwrite[1],&sa,1 << 13))
  {
    errno = EMFILE;
    return -1;
  }

  return 0;
}

#endif /* WIN32 */

/*---------------------------------------------------------------------------
  Ersatz fuer die Routine 'popen()' unter WIN32.
  ACHTUNG: Wenn 'cmd' den String '2>&1' enthaelt, wird der Standarderror File-
  Handle auf den Standardoutputfilehandle umgebogen.
  ---------------------------------------------------------------------------*/

FILE *pt_popen(const char *cmd, const char *mode)
{
#if defined(USE_CPP) && defined(_WIN32) 
  FILE *fptr = (FILE *)0;
  PROCESS_INFORMATION piProcInfo;
  STARTUPINFO siStartInfo;
  int success, umlenkung;

  my_pipein[0]   = INVALID_HANDLE_VALUE;
  my_pipein[1]   = INVALID_HANDLE_VALUE;
  my_pipeout[0]  = INVALID_HANDLE_VALUE;
  my_pipeout[1]  = INVALID_HANDLE_VALUE;
  my_pipeerr[0]  = INVALID_HANDLE_VALUE;
  my_pipeerr[1]  = INVALID_HANDLE_VALUE;

  if (mode && *mode)
  {

    my_popenmode = *mode;
    if (my_popenmode == 'r' || my_popenmode == 'w')
    {

      /*
       * Soll der stderr auf stdin umgelenkt werden ?
       */

      umlenkung = strstr("2>&1",(char *)cmd) != NULL;

      /*
       * Erzeuge die Pipes... 
       */
      if (my_pipe(my_pipein)  != -1 && my_pipe(my_pipeout) != -1  &&
         (umlenkung || my_pipe(my_pipeerr) != -1))
      {

        /*
         * Erzeuge jetzt den Sohnprozess
         */
        ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
        siStartInfo.cb           = sizeof(STARTUPINFO);
        siStartInfo.hStdInput    = my_pipein[0];
        siStartInfo.hStdOutput   = my_pipeout[1];
        siStartInfo.hStdError  = (umlenkung) ? my_pipeout[1] :  my_pipeerr[1];
        siStartInfo.dwFlags    = STARTF_USESTDHANDLES;

        success=CreateProcess(NULL,
                             (LPTSTR)cmd,       /* command line  */
                              NULL,             /* process security attributes  */
                              NULL,             /* primary thread security attributes  */
                              TRUE,             /* handles are inherited */
                              DETACHED_PROCESS, /* creation flags: Ohne Fenster (?) */
                              NULL,             /* use parent's environment */
                              NULL,             /* use parent's current directory */
                              &siStartInfo,     /* STARTUPINFO pointer */
                              &piProcInfo);     /* receives PROCESS_INFORMATION */

        if (success)
        {

          /*
           * Diese Handles gehoeren dem Sohnprozess 
           */
          CloseHandle(my_pipein[0]);  
          my_pipein[0]  = INVALID_HANDLE_VALUE;
          CloseHandle(my_pipeout[1]); 
          my_pipeout[1] = INVALID_HANDLE_VALUE;
          CloseHandle(my_pipeerr[1]); 
          my_pipeerr[1] = INVALID_HANDLE_VALUE;

          fptr = (my_popenmode=='r') ? 
                  _fdopen(_open_osfhandle((long)my_pipeout[0],_O_BINARY),"r") :
                  _fdopen(_open_osfhandle((long)my_pipein[1],_O_BINARY),"w");
        }
      }
    }
  }
  if (!fptr)
  {
    if (my_pipein[0]  != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipein[0]);
    if (my_pipein[1]  != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipein[1]);
    if (my_pipeout[0] != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipeout[0]);
    if (my_pipeout[1] != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipeout[1]);
    if (my_pipeerr[0] != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipeerr[0]);
    if (my_pipeerr[1] != INVALID_HANDLE_VALUE)
      CloseHandle(my_pipeerr[1]);
  }
  return fptr;

#elif defined(USE_CPP)  /* !WIN32  */
  return popen(cmd,mode);
#else
  return fopen(cmd,mode);
#endif /* !WIN32  */
} /* of 'pt_open' */

void pt_pclose(FILE *file)
{
#if defined(_WIN32) || !defined(USE_CPP)
  fclose(file);
#else
  pclose(file);
#endif
} /* of 'pt_close* */
