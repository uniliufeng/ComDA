/**
\file The Lorenz Model example
*/
#include "lorenz.h"
#include <iostream>
int main()
{
	ldas::Lorenz lorenz;
	double x0[3]={1.501,-1.5,25};
	for(int i=0;i<3;i++)
		std::cout<<x0[i]<<" ";
	std::cout<<std::endl;
	lorenz.parameter(x0);
	for(int i=0;i<2000;i++)
	{
		lorenz.run();
		lorenz.output();
		lorenz.parameter(lorenz.result());
	}
    return 0;
}
