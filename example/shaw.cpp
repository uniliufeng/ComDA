/**
\file The SHAW Model example
*/
#include "shaw.h"
using namespace std;

int main()
{
    ldas::SHAW shaw;
	//shaw.config('input/shaw.inp');
	shaw.run();
    return 0;
}
