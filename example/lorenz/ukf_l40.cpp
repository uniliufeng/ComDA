/**
\file UKF with Lorenz96 Model example
*/
#include "lorenz4d.h"
#include "unscentedkalmanfilter.h"
#include "random.hpp"
#include <iostream>
int main()
{
    const int steps=4000;
    const int ensemble=81;
    const int state=40;
    const int obs_delta=10;
    ldas::Lorenz4d lorenz;
    ldas::UKF ukf;
    ldas::Random random(time(NULL));
    arma::mat xtrue(state,steps),xobs(state,steps),xassim(state,steps);
    xtrue.zeros();
    xobs.zeros();
    arma::mat x0;
    x0.set_size(state,1);
    x0.fill(10);
    x0(19,0)=12;
    lorenz.parameter(x0);
    for(int i=0; i<steps; i++)
    {
        lorenz.run();

        //deal with true value
        xtrue.col(i)=lorenz.result();

        //peturb observe value
        for(int j=0; j<state; j++)
            xobs(j,i)=xtrue(j,i)+random.gaussian(0,1);
    }
    //assimilation with enkf
    ldas::Lorenz4d lz[ensemble];
    arma::mat xfen(state,ensemble),hxen(state,ensemble),yo(state,1);

    //model error, set to 0.5
    //observe error, set to 2? or 1
    arma::mat q=zeros<mat>(state,state);
    arma::mat r=zeros<mat>(state,state);
    for(int i=0; i<state; i++)
    {
        q(i,i)=1;
        r(i,i)=1;
    }
    double R=1;
    //initial state coveriance
    arma::mat P=eye<mat>(state,state);
    arma::mat x00=zeros<mat>(state,1);
    for(int i=0; i<state; i++)
    {
        x00(i,0)=x0(i,0)+random.gaussian(0,1);
    }
    arma::mat Xf;
    Xf.set_size(state,ensemble);
    //get sigma points
    ukf.sigmas(x00,P,Xf);

    for(int en=0; en<ensemble; en++)
    {
        lz[en].parameter(Xf.col(en));
    }
    for(int i=0; i<steps; i++)
    {
        for(int en=0; en<ensemble; en++)
        {
            lz[en].run();
            xfen.col(en)=lz[en].result();
            hxen.col(en)=lz[en].result();//+random.gaussian(0,1));
        }
        if((i+1)%obs_delta==0) //has observation
        {
            yo=xobs.col(i);
            ukf.update(Xf,hxen,yo,q,R,P);
        }
        else
        {
            ukf.update(Xf,q,P);
        }
        x00=ukf.getXa();
        ukf.sigmas(x00,P,Xf);

        xassim.col(i)=ukf.getXa();
        for(int en=0; en<ensemble; en++)
        {
            lz[en].parameter(Xf.col(en));
        }

        for(int j=0; j<state; j++)
            std::cout<<xtrue(j,i)<<" ";

        for(int j=0; j<state; j++)
            std::cout<<xassim(j,i)<<" ";

        double rmse=sqrt(sum(sum((xtrue.col(i)-xassim.col(i))%(xtrue.col(i)-xassim.col(i)))))/state;
        std::cout<<rmse<<" ";
        if ((i+1)%obs_delta==0)
        {
            for(int j=0; j<state; j++)
                std::cout<<xobs(j,i)<<" ";
        }
        std::cout<<std::endl;
    }
    return 0;
}
