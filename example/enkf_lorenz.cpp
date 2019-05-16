/**
\file EnKF with Lorenz Model example
*/
#include "lorenz.h"
#include "enkf.hpp"
#include "matrix.hpp"
#include <iostream>
int main()
{
    const int steps=1000;
    ldas::Lorenz lorenz;
    ldas::EnKF enkf(3,3,100);
    ldas::Random random(time(NULL));
    ldas::Matrix<double> xtrue(3,steps),xobs(3,steps),xback(3,steps),xassim(3,steps);
    xtrue.genZeros();
    xobs.genZeros();
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
        xtrue.set(1,i+1,lorenz.result()[0]);
        xtrue.set(2,i+1,lorenz.result()[1]);
        xtrue.set(3,i+1,lorenz.result()[2]);

        //peturb observe value
        xobs.set(1,i+1,lorenz.result()[0]+random.gaussian(0,1));
        xobs.set(2,i+1,lorenz.result()[1]+random.gaussian(0,1));
        xobs.set(3,i+1,lorenz.result()[2]+random.gaussian(0,1));
        lorenz.parameter(lorenz.result());
    }
    x0[0]=-5.9;
    x0[1]=-5.0;
    x0[2]=24.0;
    lorenz.parameter(x0);
    for(int i=0; i<steps; i++)
    {
        lorenz.run();
        xback.set(1,i+1,lorenz.result()[0]);
        xback.set(2,i+1,lorenz.result()[1]);
        xback.set(3,i+1,lorenz.result()[2]);
        lorenz.parameter(lorenz.result());
    }
    //assimilation with enkf
    lorenz.parameter(x0);
    ldas::Lorenz lz[100];
    ldas::Matrix<double> xfen(3,100),q(3,3),r(3,3),hxen(3,100),xaen(3,100),yo(3,1);
    double xa[3];
    //model error, set to 0.5
    //observe error, set to 2? or 1
    q.genZeros();
    r.genZeros();
    for(int i=0; i<3; i++)
    {
        q.set(i+1,i+1,1);
        r.set(i+1,i+1,1);
    }
    for(int en=0; en<100; en++)
    {
        x0[0]=-5.9+random.gaussian(0,1);
        x0[1]=-5.0+random.gaussian(0,1);
        x0[2]=24.0+random.gaussian(0,1);
        lz[en].parameter(x0);
    }
    for(int i=0; i<steps; i++)
    {
        for(int en=0; en<100; en++)
        {
            lz[en].run();
            xfen.set(1,en+1,lz[en].result()[0]);
            xfen.set(2,en+1,lz[en].result()[1]);
            xfen.set(3,en+1,lz[en].result()[2]);
            hxen.set(1,en+1,lz[en].result()[0]);//+random.gaussian(0,1));
            hxen.set(2,en+1,lz[en].result()[1]);//+random.gaussian(0,1));
            hxen.set(3,en+1,lz[en].result()[2]);//+random.gaussian(0,1));
        }
        if((i+1)%10==0) //has observation
        {
            yo=xobs.extractColV(i+1);
            enkf.update(xfen,hxen,yo,r,q);
        }
        else
        {
            enkf.update(xfen);
        }
        xaen=enkf.GetXaEn();
        xassim.setColV(i+1,enkf.GetXa());
        for(int en=0; en<100; en++)
        {
            x0[0]=xaen.get(1,en+1);
            x0[1]=xaen.get(2,en+1);
            x0[2]=xaen.get(3,en+1);
            lz[en].parameter(x0);
        }

        std::cout<<xtrue.get(1,i+1)<<" ";
        std::cout<<xtrue.get(2,i+1)<<" ";
        std::cout<<xtrue.get(3,i+1)<<" ";
        std::cout<<xback.get(1,i+1)<<" ";
        std::cout<<xback.get(2,i+1)<<" ";
        std::cout<<xback.get(3,i+1)<<" ";
        std::cout<<xassim.get(1,i+1)<<" ";
        std::cout<<xassim.get(2,i+1)<<" ";
        std::cout<<xassim.get(3,i+1)<<std::endl;
    }
    return 0;
}
