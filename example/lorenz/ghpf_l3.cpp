/**
\file GHPF with Lorenz Model example
*/
#include "lorenz.h"
#include "random.hpp"
#include "gausshermiteparticlefilter.h"
#include <iostream>
using namespace arma;
int main()
{
    const int steps=4000;
    const int ensemble=200;
    const int gpn=27;
    ldas::Lorenz lorenz;
    ldas::Random random(10);
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
    ldas::Lorenz lz[ensemble*gpn];
    ldas::GHPF ghpf(3, 3, ensemble);
    ghpf.threshold(ensemble);
    //double xa[3];
    //model error, set to 0.5
    //observe error, set to 2? or 1
    mat X_ref=zeros<mat>(3,ensemble);
    mat P[ensemble];
    for(int i=0; i<ensemble; i++)
        P[i]=eye<mat>(3,3);
    for(int en=0; en<ensemble; en++)
    {
        x0[0]=-5.9+random.gaussian(0,1);
        x0[1]=-5.0+random.gaussian(0,1);
        x0[2]=24.0+random.gaussian(0,1);
        X_ref(0,en)=x0[0];
        X_ref(1,en)=x0[1];
        X_ref(2,en)=x0[2];
    }

    mat Observation = zeros<mat>(3, 1);
    mat particles = zeros<mat>(3, ensemble*gpn);
    //mat xaen= zeros<mat>(3, ensemble);
    mat R=eye<mat>(3,3);
    mat Q=eye<mat>(3,3);
    for (int i = 0; i < steps; i++)
    {

        mat X;
        for(int en=0; en<ensemble; en++)
        {
        	//P[en]=eye<mat>(3,3);
            ghpf.gausspoint(X_ref.col(en),P[en],X);
            for(int k=0; k<gpn; k++)
            {
                x0[0]=X(0,k);
                x0[1]=X(1,k);
                x0[2]=X(2,k);
                lz[en*gpn+k].parameter(x0);
                lz[en*gpn+k].run();
                particles(0,en*gpn+k)=lz[en*gpn+k].result()[0];
                particles(1,en*gpn+k)=lz[en*gpn+k].result()[1];
                particles(2,en*gpn+k)=lz[en*gpn+k].result()[2];
            }
        }
        ghpf.update(particles,Q,P);

        if ((i+1) % 10 == 0)
        {
        	X_ref=ghpf.getXaEn();
        	//std::cout<<X_ref;
            for(int en=0; en<ensemble; en++)
            {
                ghpf.gausspoint(X_ref.col(en),P[en],X);
                for(int k=0; k<gpn; k++)
                {
                	particles.col(en*gpn+k)=X.col(k);
                }
            }
            Observation=xobs.col(i);
            /*Observation(0, 0) = xobs(0,i);
            Observation(1, 0) = xobs(1,i);
            Observation(2, 0) = xobs(2,i);*/
            //Y=X, observe equal to the state
            ghpf.update(particles,particles,Observation,R,P);
            //particle filter update
            X_ref=ghpf.getXaEn();
            //Y=X
            ghpf.update(X_ref,X_ref,Observation);
        }

        xassim.col(i)=ghpf.getXa().col(0);
        X_ref=ghpf.getXaEn();

        std::cout<<xtrue(0,i)<<" ";
        std::cout<<xtrue(1,i)<<" ";
        std::cout<<xtrue(2,i)<<" ";
        //std::cout<<xback(0,i)<<" ";
        //std::cout<<xback(1,i)<<" ";
        //std::cout<<xback(2,i)<<" ";
        std::cout<<xassim(0,i)<<" ";
        std::cout<<xassim(1,i)<<" ";
        std::cout<<xassim(2,i)<<" ";
        double rmse=sqrt(sum(sum((xtrue.col(i)-xassim.col(i))%(xtrue.col(i)-xassim.col(i)))));
        std::cout<<rmse<<" ";
        /*if ((i+1)%10==0)
        {
            std::cout<<xobs(0,i)<<" ";
            std::cout<<xobs(1,i)<<" ";
            std::cout<<xobs(2,i)<<" ";
        }*/
        std::cout<<std::endl;

    }
    return 0;
}
