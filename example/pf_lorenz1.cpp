/**
\file PF with Lorenz Model example
*/
#include "lorenz.h"
#include "random.hpp"
#include "rng.h"
#include "particlefilter.h"
#include <iostream>
int main()
{
    const int steps=1000;
    ldas::Lorenz lorenz;
    ldas::Random random(time(NULL));
    ldas::Rng rng;
    itpp::mat xtrue(3,steps),xobs(3,steps),xback(3,steps),xassim(3,steps);
    xtrue=itpp::zeros(3,steps);
    xobs=itpp::zeros(3,steps);
    double x0[3]= {-5.4458,-5.4841,22.5606};
    //for(int i=0; i<3; i++)
    //    std::cout<<x0[i]<<" ";
    //std::cout<<std::endl;
    lorenz.parameter(x0);
    for(int i=0; i<steps; i++)
    {
        lorenz.run();
        //lorenz.output();

        //deal with true value
        xtrue(0,i)=lorenz.result()[0];
        xtrue(1,i)=lorenz.result()[1];
        xtrue(2,i)=lorenz.result()[2];

        //peturb observe value
        xobs(0,i)=lorenz.result()[0]+random.gaussian(0,1);
        xobs(1,i)=lorenz.result()[1]+random.gaussian(0,1);
        xobs(2,i)=lorenz.result()[2]+random.gaussian(0,1);
        lorenz.parameter(lorenz.result());
    }
    x0[0]=-5.9;
    x0[1]=-5.0;
    x0[2]=24.0;
    lorenz.parameter(x0);
    for(int i=0; i<steps; i++)
    {
        lorenz.run();
        xback(0,i)=lorenz.result()[0];
        xback(1,i)=lorenz.result()[1];
        xback(2,i)=lorenz.result()[2];
        lorenz.parameter(lorenz.result());
    }
    //assimilation with enkf
    lorenz.parameter(x0);
    ldas::Lorenz lz[100];
    //double xa[3];
    //model error, set to 0.5
    //observe error, set to 2? or 1
    for(int en=0; en<100; en++)
    {
        x0[0]=-5.9+random.gaussian(0,1);
        x0[1]=-5.0+random.gaussian(0,1);
        x0[2]=24.0+random.gaussian(0,1);
        lz[en].parameter(x0);
    }

    itpp::mat Observation = itpp::zeros(1, 1);
    itpp::mat particles = itpp::zeros(1, 100);
    itpp::mat xaen= itpp::zeros(1, 100);
    for (int i = 0; i < steps; i++)
    {
        for(int en=0; en<100; en++)
        {
            lz[en].run();
            particles(0,en)=lz[en].result()[0];
        }


        if (i % 10 == 0)
        {
            Observation(0, 0) = xobs(0,i);
            float resampleThreshold = 3.0;
            ldas::ParticleFilter particle_filter(particles, Observation, resampleThreshold);
            particle_filter.update();
            xassim(0,i) = particle_filter.getXa()(0);
            xaen=particle_filter.getXaEn();
            for(int en=0;en<100;en++)
			{
				//每次都扰动，产生一个初值，以增加粒子的有效性
				x0[0]=xassim(0,i)*rng.uniform(0.8,1.2);
				//x0[0]=xassim(0,i);
				x0[1]=lz[en].result()[1];
				x0[2]=lz[en].result()[2];
				//x0[2]=xaen(2,en);
				lz[en].parameter(x0);
				//std::cout<<x0[0]<<" ";
			}
			//std::cout<<std::endl;
        }
        else
		{
			double x0sum=0;
			for(int en=0;en<100;en++)
			{
				x0[0]=lz[en].result()[0];
				x0[1]=lz[en].result()[1];
				x0[2]=lz[en].result()[2];
				x0sum+=x0[0];
				lz[en].parameter(x0);
				//std::cout<<x0[0]<<" ";
			}
			//std::cout<<std::endl;
			xassim(0,i) = x0sum/100.0;
		}
        double x1sum=0;
        double x2sum=0;
        for(int en=0; en<100; en++)
        {
            x1sum+=lz[en].result()[1];
            x2sum+=lz[en].result()[2];
        }
		xassim(1,i) = x1sum/100.0;
		xassim(2,i) = x2sum/100.0;
        std::cout<<xtrue(0,i)<<" ";
        std::cout<<xtrue(1,i)<<" ";
        std::cout<<xtrue(2,i)<<" ";
        std::cout<<xback(0,i)<<" ";
        std::cout<<xback(1,i)<<" ";
        std::cout<<xback(2,i)<<" ";
        std::cout<<xassim(0,i)<<" ";
        std::cout<<xassim(1,i)<<" ";
        std::cout<<xassim(2,i)<<std::endl;

    }
    return 0;
}
