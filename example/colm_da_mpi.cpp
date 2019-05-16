/**
\file The Common Land Model with EnKF data assimilation example
*/
#include "cldas_mpi.h"
using namespace std;
using namespace ldas;
int main(int argc,char** argv)
{
	Cldas cldas(argc,argv);
	cldas.config("clm.txt");
	cldas.run();
    return 0;
}
