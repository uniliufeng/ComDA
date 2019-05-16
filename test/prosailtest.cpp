#include <cpptest.h>
#include "prosail.h"
#include "radiancetransfermodel.h"

class ProsailTest : public Test::Suite
{
public:
    ProsailTest( )
    {
        TEST_ADD (ProsailTest::test1);
        TEST_ADD (ProsailTest::test2);
    }
private:
    void test1()
    {
        ldas::CanopyParameter cp;
        cp.observer_zenith=0;
        cp.solar_zenith=30;
        cp.azimuth=0;

        cp.chlorophyll_content=30;
        cp.carotenoid_content=10;
        cp.brownpigment_content=0;
        cp.cw=0.015;
        cp.cm=0.009;
        cp.structure_coefficient=1.2;
        cp.lai=2;
        cp.angle=50;
        cp.soil_coefficient=1;
        cp.skyl=70;
        cp.hot_spot=0.2;

        ldas::ProSail sail(cp);
        sail.run();
        ldas::ProSail sl(sail);
        sl.run();
        ldas::ProSail ns;
        ns.parameter(cp);
        ns.run();
        ns.parameter(sail.parameter());
        ns.run();
        //sail.output();
    }
    void test2( )
    {

    }
};

int main ( )
{
    ProsailTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

