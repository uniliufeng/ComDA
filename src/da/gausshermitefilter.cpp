#include "gausshermitefilter.h"

using namespace ldas;
GaussHermiteFilter::GaussHermiteFilter():Filter(),gauss_order(3),gauss_points_num(0)
{
    //ctor
}

GaussHermiteFilter::~GaussHermiteFilter()
{
    //dtor
}

GaussHermiteFilter::GaussHermiteFilter(const GaussHermiteFilter& other)
{
    //copy ctor
}

GaussHermiteFilter& GaussHermiteFilter::operator=(const GaussHermiteFilter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

GaussHermiteFilter::GaussHermiteFilter(const unsigned int state, const unsigned int obs):Filter(state,obs),gauss_order(3)
{
    //ctor
    init();
}

void GaussHermiteFilter::init()
{
    gauss_points_num=pow(gauss_order,state_num);
    yi=zeros<mat>(gauss_order,1);
    /*yi(0)=-1.224744871391589;
    yi(1)=0;
    yi(2)=1.224744871391589;*/
    yi(0)=-sqrt(3);
    yi(1)=0;
    yi(2)=sqrt(3);
    ai=zeros<mat>(gauss_order,1);
    /*ai(0)=0.2954089751509195;
    ai(1)=1.181635900603677;
    ai(2)=0.2954089751509196;*/
    ai(0)=1.0/6;
    ai(1)=2.0/3;
    ai(2)=1.0/6;
    qi=zeros<mat>(gauss_points_num,state_num);
    wi=zeros<mat>(gauss_points_num,state_num);
    for(int i=0; i<gauss_points_num; i++)
    {
        for(int j=state_num-1; j>=0; j--)
        {
            qi(i,j)=yi((i/int(pow(gauss_order,(state_num-j-1))))%gauss_order);
            wi(i,j)=ai((i/int(pow(gauss_order,(state_num-j-1))))%gauss_order);
        }
    }
    pwi=prod(wi,1);
}

void GaussHermiteFilter::gausspoint(const mat& Xref, const mat& P, mat& X, int flag)
{
    X=zeros<mat>(Xref.n_rows,gauss_points_num);
    for(int i=0; i<gauss_points_num; i++)
        if (flag==1)
            X.col(i)=sqrt(diagvec(P))%trans(qi.row(i))+Xref.col(0);
        else
            X.col(i)=trans(chol(P))*trans(qi.row(i))+Xref.col(0);
    /*for(int i=0;i<gauss_points_num;i++)
    	for(int j=0;j<state_num;j++)
    {
    	X(j,i)=sqrt(abs(P(j,j)))*qi(i,j)+Xref(j,0);
    }*/
}

void GaussHermiteFilter::update(const mat& X,const mat& Q, mat& P)
{
    Xa=X*pwi;
    /*Xa=zeros<mat>(state_num,1);
    for(int i=0;i<gauss_points_num;i++)
    	Xa+=X.col(i)*pwi(i);*/
    P=Q;
    for(int i=0; i<gauss_points_num; i++)
    {
        P+=(X.col(i)-Xa)*trans(X.col(i)-Xa)*pwi(i);
    }
}

void GaussHermiteFilter::update(const mat& Xt,const mat& Xf,const mat& Yt,const mat& Yo,const mat& R,mat& P)
{
    mat Zk=Yt*pwi;
    mat Pxz=zeros<mat>(Xt.n_rows,Yt.n_rows);
    mat Pzz=zeros<mat>(Yt.n_rows,Yt.n_rows);
    for(int i=0; i<gauss_points_num; i++)
    {
        Pxz+=(Xt.col(i)-Xf)*trans(Yt.col(i)-Zk)*pwi(i);
        Pzz+=(Yt.col(i)-Zk)*trans(Yt.col(i)-Zk)*pwi(i);
    }
    mat K=Pxz*inv(R+Pzz);//confirm?
    Xa=Xf+K*(Yo-Zk);
    P-=K*trans(Pxz);
}

mat GaussHermiteFilter::getXa() const
{
    return Xa;
}

double GaussHermiteFilter::rmse(const mat& xtrue) const
{
	return sqrt(accu(pow(Xa-xtrue,2))/state_num);
}
