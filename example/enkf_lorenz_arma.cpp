/**
\file EnKF with Lorenz Model example
*/
#include "lorenz.h"
#include "ensemblekalmanfilter.h"
#include "matrix.hpp"
#include <iostream>
int main()
{
    const int steps=1000;
    const int ensemble=500;
    ldas::Lorenz lorenz;
    ldas::EnKF enkf(3,3,ensemble);
    ldas::Random random(time(NULL));
    arma::mat xtrue(3,steps),xobs(3,steps),xback(3,steps),xassim(3,steps);
    xtrue.zeros();
    xobs.zeros();
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
    ldas::Lorenz lz[ensemble];
    arma::mat xfen(3,ensemble),q,r,hxen(3,ensemble),xaen(3,ensemble),yo(3,1);
    double xa[3];
    //model error, set to 0.5
    //observe error, set to 2? or 1
    q.set_size(3,3);
    q.zeros();
    r.set_size(3,3);
    r.zeros();
    for(int i=0; i<3; i++)
    {
        q(i,i)=1;
        r(i,i)=1;
    }
    for(int en=0; en<ensemble; en++)
    {
        x0[0]=-5.9+random.gaussian(0,1);
        x0[1]=-5.0+random.gaussian(0,1);
        x0[2]=24.0+random.gaussian(0,1);
        lz[en].parameter(x0);
    }
    for(int i=0; i<steps; i++)
    {
        for(int en=0; en<ensemble; en++)
        {
            lz[en].run();
            xfen(0,en)=lz[en].result()[0];
            xfen(1,en)=lz[en].result()[1];
            xfen(2,en)=lz[en].result()[2];
            hxen(0,en)=lz[en].result()[0];//+random.gaussian(0,1));
            hxen(1,en)=lz[en].result()[1];//+random.gaussian(0,1));
            hxen(2,en)=lz[en].result()[2];//+random.gaussian(0,1));
        }
        if((i+1)%10==0) //has observation
        {
            yo=xobs.col(i);
            enkf.update(xfen,hxen,yo,r,q);
        }
        else
        {
            enkf.update(xfen);
        }
        xaen=enkf.GetXaEn();
        xassim.col(i)=enkf.GetXa();
        for(int en=0; en<ensemble; en++)
        {
            x0[0]=xaen(0,en);
            x0[1]=xaen(1,en);
            x0[2]=xaen(2,en);
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
