#include <cpptest.h>
#include "landgrid.h"

class GridTest : public Test::Suite
{
public:
    GridTest( )
    {
        TEST_ADD (GridTest::test1);
        TEST_ADD (GridTest::test2);
    }
private:
    void test1()
    {
        ldas::LandGrid grid;
        grid.rows(120);
        grid.cols(160);
        grid.lat_span(0.25);
        grid.lon_span(0.25);
        grid.corner(ldas::LandPoint(49.875,72.125));
        grid.mask("/home/wlx/ldas/westchina/mask_west.900s");

        std::cout<<grid.index(159).lattitude<<" "<<grid.index(159).longitude<<std::endl;
        std::cout<<grid.index(160).lattitude<<" "<<grid.index(160).longitude<<std::endl;
        std::cout<<grid.index(19199).lattitude<<" "<<grid.index(19199).longitude<<std::endl;
        std::cout<<grid[0].lattitude<<" "<<grid[0].longitude<<std::endl;
        std::cout<<grid.mask(3,58)<<std::endl;
        TEST_ASSERT(grid[0]==grid.index(3,58));
        //TEST_ASSERT(t1>t2);
        //TEST_ASSERT(t3==t4);
    }
    void test2( ) {}
};

int main ( )
{
    GridTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}

