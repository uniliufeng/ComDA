#include <cpptest.h>
#include "permittivity.h"

class PermittivityTest : public Test::Suite
{
public:
    PermittivityTest( )
    {
        TEST_ADD (PermittivityTest::PermittivityTest1);
        TEST_ADD (PermittivityTest::PermittivityTest2);
        TEST_ADD (PermittivityTest::PermittivityTest3);
    }
private:
    void PermittivityTest1()
    {
        //default constuction test
        ldas::Permittivity r(0.915);
        std::complex<double> p=r.freewater(290);
        std::cout<<p<<std::endl;
        std::cout<<"depth:"<<r.depth(p)<<std::endl;
        r.frequency(0.0012);
        p=r.drysnow(0.3,273);
        std::cout<<p<<std::endl;
        std::cout<<"depth:"<<r.depth(p)<<std::endl;
        //TEST_ASSERT(m.get(1,1)==2);
    }
    void PermittivityTest2()
    {
        ldas::Permittivity r(5.3);
        std::cout<<r.soil(0.2,300.16,2.7,1.31,20.5,8.5)<<std::endl;
        std::cout<<r.soil(0.2,300.16,1-1.31/2.7,20.5,8.5)<<std::endl;
        std::cout<<r.soil(0.2,283.16,2.7,1.31,20.5,8.5)<<std::endl;
    }
    void PermittivityTest3()
    {
    	std::cout<<"敏感性分析"<<std::endl;
        ldas::Permittivity r(1.4);
        std::cout<<"frequency:1.4"<<std::endl;
        for(float mv=0.01;mv<0.5;mv+=0.02)
		{
			std::cout<<mv<<" "<<r.soil(mv,296.16,2.66,1.4,42.5,8.5)<<std::endl;
		}
        r.frequency(5.3);
        std::cout<<"frequency:5.3"<<std::endl;
        for(float mv=0.01;mv<0.5;mv+=0.02)
		{
			std::cout<<mv<<" "<<r.soil(mv,296.16,2.66,1.4,42.5,8.5)<<std::endl;
		}
        r.frequency(10.6);
        std::cout<<"frequency:10.6"<<std::endl;
        for(float mv=0.01;mv<0.5;mv+=0.02)
		{
			std::cout<<mv<<" "<<r.soil(mv,296.16,2.66,1.4,42.5,8.5)<<std::endl;
		}
        r.frequency(18.7);
        std::cout<<"frequency:18.7"<<std::endl;
        for(float mv=0.01;mv<0.5;mv+=0.02)
		{
			std::cout<<mv<<" "<<r.soil(mv,296.16,2.66,1.4,42.5,8.5)<<std::endl;
		}
        r.frequency(18.7);
        std::cout<<"frequency:18.7, two factor."<<std::endl;
        for(float mv=0.01;mv<0.5;mv+=0.02)
		{
		    for(float t=275;t<300;t+=0.5)
                std::cout<<mv<<" "<<t<<" "<<r.soil(mv,t,2.66,1.4,42.5,8.5)<<std::endl;
		}
    }
};

int main ( )
{
    PermittivityTest tests;
    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
}
