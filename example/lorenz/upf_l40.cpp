/**
\file UPF with Lorenz96 Model example
*/
#include "lorenz4d.h"
#include "unscentedparticlefilter.h"
#include "random.hpp"
#include <iostream>
int main()
{
    const int steps=4000;
    const int ensemble=5;
    const int state=40;
    const int obs_delta=10;
    ldas::Lorenz4d lorenz;
    ldas::UPF upf(state,state,ensemble);
    upf.threshold(ensemble);
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
    //assimilation with upf
    ldas::Lorenz4d lz[ensemble*81];
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
	//initial state coveriance
    mat P[ensemble];
	for(int i=0;i<ensemble;i++)
		P[i]=eye<mat>(state,state);
    arma::mat x00=zeros<mat>(state,1);
    mat xaen= zeros<mat>(state, ensemble);
    mat R=eye<mat>(state,state);
    for(int i=0; i<state; i++)
    {
        x00(i,0)=x0(i,0)+random.gaussian(0,1);
    }
    arma::mat Xf=zeros<mat>(state,ensemble);
    //get sigma points
    //ukf.sigmas(x00,P,Xf);

    for(int en=0; en<ensemble; en++)
    {
        lz[en].parameter(Xf.col(en));
    }
    for(int i=0; i<steps; i++)
    {
    	if ((i+1) % obs_delta == 0)
        {
            mat particles=zeros<mat>(state,ensemble*81);
            mat X;
            for(int en=0; en<ensemble; en++)
            {
            	P[en]=eye<mat>(state,state);
                upf.sigmas(Xf.col(en),P[en],X);
                for(int k=0; k<81; k++)
                {
                    lz[en*k+k].parameter(X.col(k));
                    lz[en*k+k].run();
                    particles.col(en*k+k)=lz[en*k+k].result();
                }
            }
            mat Yo=xobs.col(i);
            //Y=X, observe equal to the state
            upf.update(particles,particles,Yo,R,P);
        }
        else
        {
            mat particles=zeros<mat>(state,ensemble);
            for(int en=0; en<ensemble; en++)
            {
                lz[en].parameter(Xf.col(en));
                lz[en].run();
                particles.col(en)=lz[en].result();
            }
            upf.update(particles);
        }

       	Xf=upf.getXaEn();

        xassim.col(i)=upf.getXa();


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
