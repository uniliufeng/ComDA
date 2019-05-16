/**
\file The Common Land Model with EnKF data assimilation example
*/
#include "cldas.h"
using namespace std;
using namespace ldas;
int main()
{
	Cldas cldas;
	cldas.config("clm.txt");
	cldas.run();
    return 0;
}
