/**
\file The Lorenz Model example
*/
#include "lorenz4d.h"
#include <iostream>
int main()
{
	ldas::Lorenz4d lorenz(10,0.05);
	mat x0;
	x0.set_size(40,1);
	x0.fill(10);
	x0(19,0)=12;
	std::cout<<trans(x0);
	lorenz.parameter(x0);
	lorenz.run();
	lorenz.output();
	/*for(int i=0;i<2000;i++)
	{
		lorenz.run();
		lorenz.output();
		//lorenz.parameter(lorenz.result());
	}*/
    return 0;
}
