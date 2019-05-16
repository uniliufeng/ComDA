/**
\file The Geotop Model example
*/
#include "geotop.h"
using namespace std;
using namespace ldas;
/*----------Global variables  ------------*/
T_INIT *UV;
STRINGBIN *files;
long Nl,Nr,Nc;
double NoV;
DOUBLEMATRIX *outdata_point;
DOUBLEVECTOR *outdata_basin;
char *MISSING_FILE;
int main(int argc,char *argv[])
{
    Geotop geotop;
    geotop.config("./");
    geotop.run();
    return 0;
}
