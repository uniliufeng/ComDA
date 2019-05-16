#include <cpptest.h>
#include "qhmodel.h"
#include "radiancetransfermodel.h"
#include "permittivity.h"

class QhTest : public Test::Suite
{
public:
    QhTest( )
    {
        TEST_ADD (QhTest::test1);
        TEST_ADD (QhTest::test2);
    }
private:
    void test1()
    {
        ldas::SoilParameter sp;
        sp.cl=1.0;
        sp.theta=43.9;
        sp.frequency=5.3;
        sp.rms=0.2;
        sp.Me=0.2;
        sp.Te=300.16;
        sp.porosity=1-1.31/2.7;
        sp.sand=20.5;
        sp.clay=8.5;
        sp.We=0.2;
        sp.Op=0.05;
        sp.b1=0.62;
        sp.x=-1.08;
        sp.has_vegetation(true);
        ldas::QHModel qh(sp);
        qh.run();
        std::cout<<qh.emissivity()<<std::endl;
        std::cout<<qh.tb()<<std::endl;
    }
    void test2( )
    {
        ldas::SoilParameter sp;
        sp.cl=1.0;
        sp.theta=43.9;
        sp.frequency=5.3;
        sp.rms=0.2;
        sp.Me=0.2;
        sp.Te=300.16;
        sp.porosity=1-1.31/2.7;
        sp.sand=20.5;
        sp.clay=8.5;
        sp.We=0.2;
        //sp.Op=0.05;
        //sp.b1=0.62;
        //sp.x=-1.08;
        ldas::QHModel qh(sp);
        qh.run();
        std::cout<<qh.emissivity()<<std::endl;
        std::cout<<qh.tb()<<std::endl;
    }
};

int main ( )
{
    QhTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

