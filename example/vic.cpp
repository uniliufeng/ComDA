/**
\file The Vic-Nl(version 4.1.1) Model example
*/
#include "vic.h"
#include "global.h"
#include <iostream>
int main()
{
    ldas::Vic vic;
    vic.config("global_param.Arou");
    vic.run();
    return 0;
}
