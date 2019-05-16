/**
\file UKF with Lorenz Model example
*/
#include "lorenz.h"
#include "unscentedkalmanfilter.h"
#include "random.hpp"
//#include "matrix.hpp"
#include <iostream>
int main(int argv, char *argc[])
{
    int random_seed=0;
    int obs_delta=10;
    int qscale=1;
    switch(argv)
    {
    case 2:
        obs_delta=atoi(argc[1]);
        random_seed=time(NULL);
        break;
    case 3:
        obs_delta=atoi(argc[1]);
        random_seed=atoi(argc[2]);
        break;
    case 4:
        obs_delta=atoi(argc[1]);
        random_seed=atoi(argc[2]);
        qscale=atoi(argc[3]);
        break;
    default:
        random_seed=time(NULL);
        break;
    }

    const int steps=4000;
    ldas::Lorenz lorenz;
    ldas::UKF ukf;
    ldas::Random random(random_seed);
    arma::mat xtrue(3,steps),xobs(3,steps),xassim(3,steps);
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

    //assimilation with ukf
    const int ensemble=7;
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
    double R=1;
    R*=qscale;
    //initial state coveriance
    arma::mat P=eye<mat>(3,3);
    arma::mat x_ref;
    x_ref.set_size(3,1);
    x_ref(0,0)=-5.9;//+random.gaussian(0,1);
    x_ref(1,0)=-5.0;//+random.gaussian(0,1);
    x_ref(2,0)=24.0;//+random.gaussian(0,1);
    arma::mat Xf;
    Xf.set_size(3,ensemble);
    //get sigma points
    ukf.sigmas(x_ref,P,Xf);
    for(int en=0; en<ensemble; en++)
    {
        x0[0]=Xf(0,en);
        x0[1]=Xf(1,en);
        x0[2]=Xf(2,en);
        lz[en].parameter(x0);
    }
    for(int i=0; i<steps; i++)
    {
        for(int en=0; en<ensemble; en++)
        {
            lz[en].run();
            Xf(0,en)=lz[en].result()[0];
            Xf(1,en)=lz[en].result()[1];
            Xf(2,en)=lz[en].result()[2];
            hxen(0,en)=lz[en].result()[0];//+random.gaussian(0,1));
            hxen(1,en)=lz[en].result()[1];//+random.gaussian(0,1));
            hxen(2,en)=lz[en].result()[2];//+random.gaussian(0,1));
        }
        if((i+1)%obs_delta==0) //has observation
        {
            mat yo=xobs.col(i);
            ukf.update(Xf,hxen,yo,q,R,P);
        }
        else
        {
            ukf.update(Xf,q,P);
        }
        x_ref=ukf.getXa();
        ukf.sigmas(x_ref,P,Xf);
        for(int en=0; en<ensemble; en++)
        {
            x0[0]=Xf(0,en);
            x0[1]=Xf(1,en);
            x0[2]=Xf(2,en);
            lz[en].parameter(x0);
        }

        std::cout<<xtrue(0,i)<<" ";
        std::cout<<xtrue(1,i)<<" ";
        std::cout<<xtrue(2,i)<<" ";
        std::cout<<x_ref(0,0)<<" ";
        std::cout<<x_ref(1,0)<<" ";
        std::cout<<x_ref(2,0)<<" ";
        std::cout<<ukf.rmse(xtrue.col(i))<<" ";
        if ((i+1)%obs_delta==0)
        {
            std::cout<<xobs(0,i)<<" ";
            std::cout<<xobs(1,i)<<" ";
            std::cout<<xobs(2,i)<<" ";
            std::cout<<ukf.rmse(xobs.col(i))<<" ";
        }
        std::cout<<std::endl;
    }
    return 0;
}
