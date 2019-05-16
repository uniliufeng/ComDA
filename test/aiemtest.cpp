#include <cpptest.h>
#include "aiem.h"
#include "radiancetransfermodel.h"
#include "permittivity.h"

class AiemTest : public Test::Suite
{
public:
    AiemTest( )
    {
        TEST_ADD (AiemTest::test1);
        TEST_ADD (AiemTest::test2);
    }
private:
    void test1()
    {
        ldas::SoilParameter sp;
        sp.cl=3;
        sp.theta=50;
        sp.frequency=19;
        sp.rms=1;
        ldas::AIEM aiem(sp);
        aiem.permittivity(std::complex<double>(5.1,0.4));
        //aiem.gauss_integral(16,16);
        aiem.run();
        ldas::Emissivity e=aiem.emissivity();
        std::cout<<e<<std::endl;
        //TEST_ASSERT(grid[0]==grid.index(3,58));
        //TEST_ASSERT(t1>t2);
        //TEST_ASSERT(t3==t4);
    }
    void test2( )
    {
        ldas::SoilParameter sp;
        sp.cl=1.0;
        sp.theta=43.9;
        sp.frequency=5.3;
        sp.rms=0.2;
        ldas::AIEM aiem(sp);
        //aiem.permittivity(p.soil(0.4,300.16,1-1.31/2.7,0.205,0.085));
        aiem.permittivity(std::complex<double>(7.5968899726867676,0.81962758302688599));
        aiem.fresnel(ldas::AIEM_SIMULTANEOUS);
        aiem.corelation(ldas::AIEM_EXPONENT);
        aiem.mode(ldas::AIEM_ACTIVE);
        //aiem.gauss_integral(16,16);
        aiem.run();
        std::cout<<aiem.scattering()<<std::endl;
    }
};

int main ( )
{
    AiemTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

