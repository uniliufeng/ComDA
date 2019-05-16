#include <cpptest.h>
#include "memls.h"
#include "radiancetransfermodel.h"
#include "matrix.hpp"

class MemlsTest : public Test::Suite
{
public:
    MemlsTest( )
    {
        TEST_ADD (MemlsTest::test1);
        TEST_ADD (MemlsTest::test2);
    }
private:
    void test1()
    {
        ldas::SnowParameter sp;
        sp.theta=10;
        sp.frequency=10.87;
        sp.col=1;
        sp.temperature_sky=0;
        sp.temperature_ground=273;
        sp.depth.setDim(1,sp.col);
        sp.depth.set(1,sp.col,66);
        sp.density.setDim(1,sp.col);
        sp.density.set(1,sp.col,0.3);
        sp.temperature.setDim(1,sp.col);
        sp.temperature.set(1,sp.col,273);
        sp.wetness.setDim(1,sp.col);
        sp.wetness.set(1,sp.col,0.06);
        sp.pci.setDim(1,sp.col);
        sp.pci.set(1,sp.col,0.12);
        ldas::MEMLS memls(sp);
        memls.run();
        std::cout<<memls.brightness()<<std::endl;
        //TEST_ASSERT(grid[0]==grid.index(3,58));
        //TEST_ASSERT(t1>t2);
        //TEST_ASSERT(t3==t4);
    }
    void test2( )
    {
    }
};

int main ( )
{
    MemlsTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

