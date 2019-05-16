/**
\file PF with Lorenz Model example
*/
#include "lorenz.h"
#include "random.hpp"
#include "particlefilter.h"
#include <iostream>
using namespace arma;
int main()
{
    const int steps=1000;
    const int ensemble=500;
    ldas::Lorenz lorenz;
    ldas::Random random(time(NULL));
    mat xtrue=zeros<mat>(3,steps);
    mat xobs=zeros<mat>(3,steps);
    mat xback=zeros<mat>(3,steps);
    mat xassim=zeros<mat>(3,steps);
    double x0[3]= {-5.4458,-5.4841,22.5606};
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
    ldas::Lorenz lz[ensemble];
    //double xa[3];
    //model error, set to 0.5
    //observe error, set to 2? or 1
    for(int en=0; en<ensemble; en++)
    {
        x0[0]=-5.9+random.gaussian(0,1);
        x0[1]=-5.0+random.gaussian(0,1);
        x0[2]=24.0+random.gaussian(0,1);
        lz[en].parameter(x0);
    }

    mat Observation = zeros<mat>(3, 1);
    mat particles = zeros<mat>(3, ensemble);
    mat xaen= zeros<mat>(3, ensemble);
    ldas::ParticleFilter particle_filter(3, 3, ensemble);
    particle_filter.threshold(1.0*ensemble);
    //particle_filter.threshold(5.0);
    for (int i = 0; i < steps; i++)
    {
        for(int en=0; en<ensemble; en++)
        {
            lz[en].run();
            particles(0,en)=lz[en].result()[0];
            particles(1,en)=lz[en].result()[1];
            particles(2,en)=lz[en].result()[2];
        }

        if ((i+1) % 10 == 0)
        {
            Observation(0, 0) = xobs(0,i);
            Observation(1, 0) = xobs(1,i);
            Observation(2, 0) = xobs(2,i);
            // Yf=Xf here
            particle_filter.update(particles,particles,Observation);
        }
        else
        {
            particle_filter.update(particles);
        }
        xassim(0,i) = particle_filter.getXa()(0);
        xassim(1,i) = particle_filter.getXa()(1);
        xassim(2,i) = particle_filter.getXa()(2);
        for(int en=0; en<ensemble; en++)
        {
            //每次都扰动，产生一个初值，以增加粒子的有效性
            if ((i+1) % 100 == 0)
            {
                x0[0]=xassim(0,i)*random.uniform(0.8,1.2);
                x0[1]=xassim(1,i)*random.uniform(0.8,1.2);
                x0[2]=xassim(2,i)*random.uniform(0.8,1.2);
            }
            else
            {
                x0[0]=particle_filter.getXaEn()(0,en);//lz[en].result()[0];
                x0[1]=particle_filter.getXaEn()(1,en);//lz[en].result()[1];
                x0[2]=particle_filter.getXaEn()(2,en);//lz[en].result()[2];
            }
            lz[en].parameter(x0);
        }

        std::cout<<xtrue(0,i)<<" ";
        std::cout<<xtrue(1,i)<<" ";
        std::cout<<xtrue(2,i)<<" ";
        std::cout<<xback(0,i)<<" ";
        std::cout<<xback(1,i)<<" ";
        std::cout<<xback(2,i)<<" ";
        std::cout<<xassim(0,i)<<" ";
        std::cout<<xassim(1,i)<<" ";
        std::cout<<xassim(2,i)<<" ";
        double rmse=sqrt(sum(sum((xtrue.col(i)-xassim.col(i))%(xtrue.col(i)-xassim.col(i)))));
        std::cout<<rmse<<" ";
        if ((i+1)%10==0)
        {
            std::cout<<xobs(0,i)<<" ";
            std::cout<<xobs(1,i)<<" ";
            std::cout<<xobs(2,i)<<" ";
        }
        std::cout<<std::endl;

    }
    return 0;
}
