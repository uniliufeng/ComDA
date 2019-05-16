/**
\file The Noah LSM Model example
*/
#include "noahlsm.h"
using namespace std;
using namespace ldas;

int main()
{
    Noahlsm noah;
    noah.config("bondville.dat");
    noah.run();
    return 0;
}
