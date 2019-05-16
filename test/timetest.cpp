#include <cpptest.h>
#include "time.hpp"

class TimeTest : public Test::Suite
{
public:
    TimeTest()
    {
        TEST_ADD(TimeTest::test1);
        TEST_ADD(TimeTest::test2);
        TEST_ADD(TimeTest::test3);
        TEST_ADD(TimeTest::test4);
    }
private:
	//test time addXXX methods.
    void test1()
    {
        ldas::Time t1(2010,3,28,10);
        t1.addDays(4);
        t1.addHours(4);
        t1.addMinutes(4);
        t1.addSeconds(2000);
        ldas::Time t2(2010,4,1,14,37,20);
        TEST_ASSERT(t1==t2);
        t1.addSeconds(-2000);
        t1.addMinutes(-4);
        t1.addHours(-4);
        t1.addDays(-4);
        ldas::Time t3(2010,3,28,10);
        TEST_ASSERT(t1==t3);
    }
    //test addDays with big number
    void test2()
    {
        ldas::Time t(2010,1,1);
        t.addSeconds(-10);
        ldas::Time t1(2009,12,31,23,59,50);
        TEST_ASSERT(t==t1);
        t.julian(2010,1);
        t.addDays(-1462);
        t1.julian(2005,365);
        TEST_ASSERT(t==t1);
        t.addDays(1462);
        t1.julian(2010,1);
        TEST_ASSERT(t==t1);
        t.addDays(50);
        t.addDays(-50);
        TEST_ASSERT(t==t1);
    }
    //test ISO8601 construction
    void test3()
    {
        ldas::Time t1("1998-01-01 06:30");
        ldas::Time t2(1998,1,1,6,30);
        TEST_ASSERT(t1==t2);
        ldas::Time t3("1997-07-16T19:20:30+01:00");
        ldas::Time t4(1997,07,16,19,20,30,1);
        TEST_ASSERT(t3==t4);
    }
    //test different years
    void test4()
    {
        ldas::Time t1(2007,12,31,23);
        t1.addHours(1);
        ldas::Time t2(2008,1,1);
        TEST_ASSERT(t1==t2);
        ldas::Time t3(2007,12,31);
        t3.addDays(1);
        TEST_ASSERT(t1==t3);
        ldas::Time t4(2008,12,31);
        t4.addDays(1);
        std::cout<<t4;
    }
};

int main ( )
{
    TimeTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}
