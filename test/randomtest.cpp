#include <cpptest.h>
#include "random.hpp"

class RandomTest : public Test::Suite
{
public:
    RandomTest( )
    {
        TEST_ADD (RandomTest::RandomTest1);
        TEST_ADD (RandomTest::RandomTest2);
        TEST_ADD (RandomTest::RandomTest3);
    }
private:
    void RandomTest1()
    {
        //default constuction test
        ldas::Random r;
        std::cout<<r.gaussian(0,1)<<" "<<r.gaussian(0,1)<<std::endl;
        //TEST_ASSERT(m.get(1,1)==2);
    }
    void RandomTest2()
    {

    }
    void RandomTest3()
    {

    }
};

int main ( )
{
    RandomTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}
