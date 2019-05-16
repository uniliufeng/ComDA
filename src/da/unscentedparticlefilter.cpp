#include "unscentedparticlefilter.h"
using namespace ldas;

UnscentedParticleFilter::UnscentedParticleFilter():m_alpha(1),m_beta(2),m_kappa(0),ParticleFilter()
{
    //ctor
}

UnscentedParticleFilter::~UnscentedParticleFilter()
{
    //dtor
}

UnscentedParticleFilter::UnscentedParticleFilter(const UnscentedParticleFilter& other)
{
    //copy ctor
}

UnscentedParticleFilter& UnscentedParticleFilter::operator=(const UnscentedParticleFilter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

UnscentedParticleFilter::UnscentedParticleFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum, const double a, const double b, const double k)
    :ParticleFilter(StateNum,ObsNum,EnNum),m_alpha(a),m_beta(b),m_kappa(k)
{
    //ctor
}

void UnscentedParticleFilter::update(const mat& state_en/*_sigmas, mat* P*/)
{
    /*int sigma_num=state_en_sigmas.n_rows*2+1;
    mat state_en=zeros<mat>(state_en_sigmas.n_rows,particle_num);
    for(int k=0; k<particle_num; k++)
    {
        //sigma points,get X_ref, computed in the foreign program
        //model computation,get Xf
        //ukf update, get Xa
        mat Xf=zeros<mat>(state_en_sigmas.n_rows,sigma_num);
        for(int i=0; i<sigma_num; i++)
        {
            Xf.col(i)=state_en_sigmas.col(k*sigma_num+i);
        }
        mat Xa_p=zeros<mat>(state_en_sigmas.n_rows,1);
        ukf(Xf,Xa_p,P[k]);
        state_en.col(k)=Xa_p.col(0);
    }*/
    particles=state_en;
    //Xa=mean(particles,1);
    mat weights_en=zeros<mat>(state_num,particle_num);
    for(int i=0; i<state_num; i++)
        weights_en.row(i)=weights.row(0);
    Xa=sum(weights_en%particles,1);
}

void UnscentedParticleFilter::update(const mat& state_en_sigmas, const mat& obs_forecast_en_sigmas, const mat& obs, const mat& R, mat* P)
{
    int sigma_num=state_en_sigmas.n_rows*2+1;
    mat state_en=zeros<mat>(state_en_sigmas.n_rows,particle_num);
    for(int k=0; k<particle_num; k++)
    {
        //sigma points,get X_ref, computed in the foreign program
        //model computation,get Xf
        //ukf update, get Xa
        mat Xf=zeros<mat>(state_en_sigmas.n_rows,sigma_num);
        mat Yf=zeros<mat>(obs_forecast_en_sigmas.n_rows,sigma_num);
        for(int i=0; i<sigma_num; i++)
        {
            Xf.col(i)=state_en_sigmas.col(k*sigma_num+i);
            Yf.col(i)=obs_forecast_en_sigmas.col(k*sigma_num+i);
        }
        mat Xa_p=zeros<mat>(state_en_sigmas.n_rows,1);
        ukf(Xf,Yf,obs,R,Xa_p,P[k]);
        state_en.col(k)=Xa_p.col(0);
    }
    mat obs_en=zeros<mat>(observe_num,particle_num);
    for(int i=0; i<particle_num; i++)
        obs_en.col(i)=obs.col(0);
    mat llh=likelihood(obs_en,state_en)+1e-99;
    weights=weights%llh;
    weights=weights/sum(sum(weights));

    // calculate effective particle set size
    double Neff = 1.0 / sum(sum(pow(weights, 2)));
    if (Neff < m_threshold)
    {
        particles=residual_resample(weights,state_en);
        weights = 1.0/particle_num*ones<mat>(1,particle_num);
    }
    else
    {
        particles=state_en;
    }
    mat weights_en=zeros<mat>(state_num,particle_num);
    for(int i=0; i<state_num; i++)
        weights_en.row(i)=weights.row(0);
    Xa=sum(weights_en%particles,1);
}
void UnscentedParticleFilter::ukf(const mat& Xf, const mat& Yf, const mat& Yo, const mat& R, mat& Xa, mat& P)
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
    //mat Q=eye(Xf.n_rows,Xf.n_rows);
    P1=X2*diagmat(Wc.row(0))*trans(X2);//+Q;
    //unscented transformation of measurments
    mat z1=zeros<mat>(observe_num,1);
    mat Z2=zeros<mat>(observe_num,Yf.n_cols);
    mat P2=zeros<mat>(observe_num,observe_num);
    ut(Yf,Wm,Wc,z1,Z2);
    //mat R=eye(Yf.n_rows,Yf.n_rows);
    P2=Z2*diagmat(Wc.row(0))*trans(Z2)+R;
    //if (det(P2)<1e-10) P2=P2+R;

    //transformed cross-covariance
    mat P12=X2*diagmat(Wc.row(0))*trans(Z2);
    mat K=P12*inv(P2);
    //state update
    Xa=x1+K*(Yo-z1);
    //covariance update
    P=P1-K*trans(P12);
}

void UnscentedParticleFilter::ukf(const mat& Xf, mat& X, mat& P)
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
    ut(Xf,Wm,Wc,X,X2);
    //mat Q=eye<mat>(Xf.n_rows,Xf.n_rows);
    P=X2*diagmat(Wc.row(0))*trans(X2);//+Q;
}
void UnscentedParticleFilter::ut(const mat& Y,const mat& Wm,const mat& Wc, mat& y, mat& Y1)
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
void UnscentedParticleFilter::sigmas(const mat& x,const mat& P,  mat& X)
{
    state_num=x.n_rows;
    /** c: coefficient */
    double c=m_alpha*m_alpha*(state_num+m_kappa);
    c=sqrt(c);
    mat A = c*trans(chol(P));
    //copy x.n_elem columns
    mat Y;
    Y.set_size(x.n_rows,x.n_elem);
    for(int i=0; i<x.n_elem; i++)
        Y.col(i)=x.col(0);
    X=join_rows(x,Y+A);
    X=join_rows(X,Y-A);
}
