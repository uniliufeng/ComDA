#include "unscentedkalmanfilter.h"

using namespace ldas;
UnscentedKalmanFilter::UnscentedKalmanFilter():Filter(),m_alpha(1.0),m_beta(2.0),m_kappa(0)
{
    //ctor
}

UnscentedKalmanFilter::~UnscentedKalmanFilter()
{
    //dtor
}

UnscentedKalmanFilter::UnscentedKalmanFilter(const UnscentedKalmanFilter& other)
{
    //copy ctor
}

UnscentedKalmanFilter& UnscentedKalmanFilter::operator=(const UnscentedKalmanFilter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

UnscentedKalmanFilter::UnscentedKalmanFilter(const unsigned int state, const unsigned int obs, const double a, const double b, const double k)
    :Filter(state,obs),m_alpha(a),m_beta(b),m_kappa(k)
{
    //ctor
}
void UnscentedKalmanFilter::update(const mat& Xf,const mat& Yf, const mat& Yo, const mat& Q, const double R, mat& P)
{
    state_num=Xf.n_rows;
    observe_num=Yf.n_rows;
    double lambda=m_alpha*m_alpha*(state_num+m_kappa)-state_num;
    double c=state_num+lambda;
    mat Wm=zeros<mat>(1,state_num*2+1);
    Wm(0)=lambda/c;
    for(int i=1; i<=state_num*2; i++)
        Wm(i)=0.5/c;
    mat Wc=Wm;
    Wc(0)=Wc(0)+(1-m_alpha*m_alpha+m_beta);
    //unscented transformation of process
    mat x1=zeros<mat>(state_num,1);
    mat X2=zeros<mat>(state_num,Xf.n_cols);
    mat P1=zeros<mat>(state_num,state_num);
    ut(Xf,Wm,Wc,x1,X2);
    P1=X2*diagmat(Wc.row(0))*trans(X2)+Q;
    //unscented transformation of measurments
    mat z1=zeros<mat>(observe_num,1);
    mat Z2=zeros<mat>(observe_num,Yf.n_cols);
    mat P2=zeros<mat>(observe_num,observe_num);
    ut(Yf,Wm,Wc,z1,Z2);
    P2=Z2*diagmat(Wc.row(0))*trans(Z2)+R;
    //transformed cross-covariance
    mat P12=X2*diagmat(Wc.row(0))*trans(Z2);
    mat K=P12*inv(P2);
    //state update
    Xa=x1+K*(Yo-z1);
    //covariance update
    P=P1-K*trans(P12);
}
void UnscentedKalmanFilter::update(const mat& Xf, const mat& Q, mat& P)
{
    state_num=Xf.n_rows;
    Xa.set_size(state_num,1);
    Xa.zeros();
    double lambda=m_alpha*m_alpha*(state_num+m_kappa)-state_num;
    double c=state_num+lambda;
    mat Wm=zeros<mat>(1,state_num*2+1);
    Wm(0)=lambda/c;
    for(int i=1; i<=state_num*2; i++)
        Wm(i)=0.5/c;
    mat Wc=Wm;
    Wc(0)=Wc(0)+(1-m_alpha*m_alpha+m_beta);
    //unscented transformation of process
    //mat x1=zeros<mat>(state_num,1);
    mat X2=zeros<mat>(state_num,Xf.n_cols);
    //mat P1=zeros<mat>(state_num,state_num);
    ut(Xf,Wm,Wc,Xa,X2);
    P=X2*diagmat(Wc.row(0))*trans(X2)+Q;
}

/** Unscented Transformation */
void UnscentedKalmanFilter::ut(const mat& Y,const mat& Wm,const mat& Wc, mat& y, mat& Y1)
{
    unsigned int L=Y.n_cols;
    //y=zeros<mat>(n,1);
    //mat Y=zeros<mat>(n,L);
    for(int k=1; k<=L; k++)
    {
        //Y.col(k-1)=f(X.col(k-1));
        y=y+Wm(k-1)*Y.col(k-1);
    }
    Y1=zeros<mat>(Y.n_rows,L);
    for(int i=0; i<L; i++)
        Y1.col(i)=Y.col(i)-y.col(0);
    //P=Y1*diag(Wc)*htrans(Y1)+R;
}

/**  Sigma points around reference point */
void UnscentedKalmanFilter::sigmas(const mat& x,const mat& P,  mat& X)
{
    state_num=x.n_rows;
    /** c: coefficient */
    double c=m_alpha*m_alpha*(state_num+m_kappa);
    c=sqrt(c);
    mat A = c*trans(chol(P));
    //复制x.n_elem列
    mat Y;
    Y.set_size(x.n_rows,x.n_elem);
    for(int i=0; i<x.n_elem; i++)
        Y.col(i)=x.col(0);
    X=join_rows(x,Y+A);
    X=join_rows(X,Y-A);
}

mat UnscentedKalmanFilter::getXa() const
{
    return Xa;
}

double UnscentedKalmanFilter::rmse(const mat& xtrue) const
{
	return sqrt(accu(pow(Xa-xtrue,2))/state_num);
}
