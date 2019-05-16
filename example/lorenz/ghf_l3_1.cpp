/**
\file GHF with Lorenz Model example, one observation
*/
#include "lorenz.h"
#include "random.hpp"
#include "gausshermitefilter.h"
#include <iostream>
#include <cstdlib>

using namespace arma;
int main(int argv, char *argc[])
{
	int random_seed=0;
	int obs_delta=10;
	int use_sqrt=1;
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
	default:
		random_seed=time(NULL);
		break;
	}

    const int steps=4000;
    const int ensemble=27;
    const int observe_num=2;
    ldas::Lorenz lorenz;
    ldas::Random random(random_seed);
    mat xtrue=zeros<mat>(3,steps);
    mat xobs=zeros<mat>(3,steps);
    mat xassim=zeros<mat>(3,steps);
    double x0[3]= {-5.4458,-5.4841,22.5606};
    //double x0[3]= {1.508870,-1.531271,25.46091};
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

    ldas::Lorenz lz[ensemble];

    ldas::GHF ghf(3, observe_num);
    mat Xt=zeros<mat>(3,ensemble);
    mat P=eye<mat>(3,3);
    mat Q=eye<mat>(3,3);
    /*Q(0,0)=1.0;
    Q(1,1)=1.0;
    Q(2,2)=15.0;*/
    P=Q;
    mat R=eye<mat>(observe_num,observe_num);
    //R=R*3;
    mat X_ref=zeros<mat>(3,1);
    X_ref(0)=-5.9;//+random.gaussian(0,1);
    X_ref(1)=-5;//+random.gaussian(0,1);
    X_ref(2)=24;//+random.gaussian(0,1);

    for (int i = 0; i < steps; i++)
    {
    	ghf.gausspoint(X_ref,P,Xt,use_sqrt);
        for(int en=0; en<ensemble; en++)
        {
            x0[0]=Xt(0,en);
            x0[1]=Xt(1,en);
            x0[2]=Xt(2,en);
            lz[en].parameter(x0);
            lz[en].run();
            Xt(0,en)=lz[en].result()[0];
            Xt(1,en)=lz[en].result()[1];
            Xt(2,en)=lz[en].result()[2];
        }
        ghf.update(Xt,Q,P);

        if ((i+1) % obs_delta == 0)
        {
            X_ref=ghf.getXa();
            //std::cout<<X_ref;
            ghf.gausspoint(X_ref,P,Xt,use_sqrt);
            mat Yo=zeros<mat>(observe_num,1);
            Yo(0)=xobs(0,i);
            Yo(1)=xobs(1,i);
            mat Yt=zeros<mat>(observe_num,ensemble);
            Yt.row(0)=Xt.row(0);
            Yt.row(1)=Xt.row(1);
            // Yt=Xt here
            ghf.update(Xt,X_ref,Yt,Yo,R,P);
        }
        X_ref=ghf.getXa();
        xassim.col(i) = X_ref;

        std::cout<<xtrue(0,i)<<" ";
        std::cout<<xtrue(1,i)<<" ";
        std::cout<<xtrue(2,i)<<" ";

        std::cout<<xassim(0,i)<<" ";
        std::cout<<xassim(1,i)<<" ";
        std::cout<<xassim(2,i)<<" ";
        //double rmse=sqrt(sum(sum((xtrue.col(i)-xassim.col(i))%(xtrue.col(i)-xassim.col(i)))));
        //double rmse=sqrt(pow(xtrue(0,i)-xassim(0,i),2));
        std::cout<<ghf.rmse(xtrue.col(i))<<" ";
        if ((i+1)%obs_delta==0)
        {
            std::cout<<xobs(0,i)<<" ";
            std::cout<<xobs(1,i)<<" ";
            std::cout<<xobs(2,i)<<" ";
            double rmsd=sqrt((pow(xassim(0,i)-xobs(0,i),2)+pow(xassim(1,i)-xobs(1,i),2))/observe_num);
            std::cout<<rmsd<<" ";
        }
        std::cout<<std::endl;

    }
    //compute NSE
    mat avg_q=zeros<mat>(observe_num,1);
    mat nse1=zeros<mat>(observe_num,1);
    mat nse2=zeros<mat>(observe_num,1);
    int count=0;
    for(int i=0;i<steps;i++)
	{
		if ((i+1) % obs_delta == 0)
		{
			avg_q(0)+=xobs(0,i);
			avg_q(1)+=xobs(1,i);
			count++;
		}
	}
	avg_q/=count;
	for(int i=0;i<steps;i++)
	{
		if ((i+1) % obs_delta == 0)
		{
			nse1(0)+=pow(xassim(0,i)-xobs(0,i),2);
			nse1(1)+=pow(xassim(1,i)-xobs(1,i),2);
			nse2(0)+=pow(xobs(0,i)-avg_q(0),2);
			nse2(1)+=pow(xobs(1,i)-avg_q(1),2);
		}
	}
	mat nse=1-nse1/nse2;
	std::cout<<"#NSE "<<nse(0)<<" "<<nse(1)<<std::endl;
    return 0;
}
