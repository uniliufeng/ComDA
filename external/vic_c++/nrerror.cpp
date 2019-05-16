#include "vic.h"
using namespace ldas;

static char vcid[] = "$Id: nrerror.c,v 3.2 1999/08/23 23:59:06 vicadmin Exp $";

void Vic::nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	//void _exit();

	fprintf(stderr,"Model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
