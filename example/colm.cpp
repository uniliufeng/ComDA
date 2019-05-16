/**
\file The Common Land Model example
*/
#include "commonlandmodel.h"
using namespace std;
using namespace ldas;
using namespace ldas::colm;
int main()
{
    CommonLandModel colm;
    try
    {
        colm.config("clm.txt");
        colm.run();
    }
    catch(Exception e)
    {
        e.what();
    }
    return 0;
}
