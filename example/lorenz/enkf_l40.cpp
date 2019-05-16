/**
\file EnKF with Lorenz96 Model example
*/
#include "lorenz4d.h"
#include "ensemblekalmanfilter.h"
#include <iostream>
int main()
{
    const int steps=4000;
    const int ensemble=400;
    const int state=40;
    const int obs_delta=10;
    ldas::Lorenz4d lorenz;
    ldas::EnKF enkf(state,state,ensemble);
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
        for(int j=0;j<state;j++)
			xobs(j,i)=xtrue(j,i)+random.gaussian(0,1);
    }
    //assimilation with enkf
    ldas::Lorenz4d lz[ensemble];
    arma::mat xfen(state,ensemble),hxen(state,ensemble),xaen(state,ensemble),yo(state,1);

    //model error, set to 0.5
    //observe error, set to 2? or 1
    arma::mat q=zeros<mat>(state,state);
    arma::mat r=zeros<mat>(state,state);
    for(int i=0; i<state; i++)
    {
        q(i,i)=1;
        r(i,i)=1;
    }
    arma::mat x00=zeros<mat>(state,1);
    for(int en=0; en<ensemble; en++)
    {
    	for(int i=0;i<state;i++)
		{
			x00(i,0)=x0(i,0)+random.gaussian(0,1);
		}
        lz[en].parameter(x00);
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
            lz[en].parameter(xaen.col(en));
        }

        for(int j=0;j<state;j++)
			std::cout<<xtrue(j,i)<<" ";

        for(int j=0;j<state;j++)
			std::cout<<xassim(j,i)<<" ";

        double rmse=sqrt(sum(sum((xtrue.col(i)-xassim.col(i))%(xtrue.col(i)-xassim.col(i)))))/state;
        std::cout<<rmse<<" ";
        if ((i+1)%obs_delta==0)
        {
        	for(int j=0;j<state;j++)
				std::cout<<xobs(j,i)<<" ";
        }
        std::cout<<std::endl;
    }
    return 0;
}
