#include "ensemblekalmanfilter.h"

using namespace ldas;
EnsembleKalmanFilter::EnsembleKalmanFilter():Filter(),ensemble_num(0),m_perturb_model(true),m_perturb_observe(true)
{
    //ctor
}

EnsembleKalmanFilter::~EnsembleKalmanFilter()
{
    //dtor
}

EnsembleKalmanFilter::EnsembleKalmanFilter(const EnsembleKalmanFilter& other)
{
    //copy ctor
}

EnsembleKalmanFilter& EnsembleKalmanFilter::operator=(const EnsembleKalmanFilter& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

EnsembleKalmanFilter::EnsembleKalmanFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum):ensemble_num(EnNum),m_perturb_model(true),m_perturb_observe(true)
{
    //ctor
    state_num=StateNum;
    observe_num=ObsNum;
    Xf.set_size(state_num, 1);
    Xa.set_size(state_num, 1);
    XfEn.set_size(state_num, ensemble_num);
    XaEn.set_size(state_num, ensemble_num);
    //Y.set_size(observe_num, 1);
    YEn.set_size(observe_num, ensemble_num);
    Yo.set_size(observe_num, 1);
    YoEn.set_size(observe_num, ensemble_num);
    //R.set_size(observe_num, observe_num);
    Pf.set_size(state_num, state_num);
    Pa.set_size(state_num, state_num);
    K.set_size(state_num, observe_num);
}

void EnsembleKalmanFilter::update(mat& m_XfEn, mat& m_YEn,
                                  mat& m_Yo, mat& m_R, mat& m_Q)
{
    state_num  = m_XfEn.n_rows;
    ensemble_num  = m_XfEn.n_cols;
    observe_num  = m_Yo.n_rows;

    //Q=m_Q;
    XfEn = m_XfEn;
    YEn  = m_YEn;
    Yo   = m_Yo;
    R    = m_R;

    YoEn.set_size(observe_num,ensemble_num);
    Xa.set_size(state_num,1);
    Xf.set_size(state_num,1);
    XaEn.set_size(state_num,ensemble_num);
    Pa.set_size(state_num,state_num);
    Pf.set_size(state_num,state_num);

    //mat	PH(state_num, observe_num);
    //mat  HPH(observe_num, observe_num);
    //mat	XfEnInc(state_num, ensemble_num);
    //mat	YEnInc(observe_num, ensemble_num);
    //mat	D(observe_num, ensemble_num);			// ensemble of inovation

    //mat  temp(observe_num,state_num);

    //generate XfEn with perturb model error
    if (m_perturb_model)
        predict(m_Q);

    //perturb observation
    if (m_perturb_observe)
        perturbYo(m_R);
    else
        for(int i=0; i<ensemble_num; i++)
            YoEn.col(i)=Yo;
    //D= YoEn - YEn;

    // increment
    // Ensemble perturbation matrix of model forecast
    mat XfEnInc	= ensembleInc(XfEn);
    mat YEnInc	= ensembleInc(YEn);
    //mat tmp(ensemble_num,observe_num);
    //tmp=trans(YEnInc);

    // Pf*H()
    mat PH		= XfEnInc*trans(YEnInc)/(ensemble_num-1.0);
    // H()*Pf*H()^T
    mat HPH		= YEnInc *trans(YEnInc)/(ensemble_num-1.0);
    //temp = HPH+R;

    try
    {
        K		= PH*inv(HPH+R);
        XaEn	= XfEn + K*(YoEn-YEn);
    }
    catch (...)
    {
        XaEn	= XfEn;
    }
    if (Xl.n_rows==state_num && Xh.n_rows==state_num)
        for (int i=0; i<ensemble_num; i++)
        {
            mat Xai	= XaEn.col(i);
            for (int j=0; j<state_num; j++)
            {
                if ( Xai(j, 0) < Xl(j, 0) )
                    Xai(j, 0)=Xl(j,0);
                if ( Xai(j, 0) > Xh(j, 0) )
                    Xai(j, 0)=Xh(j, 0);
            }
            XaEn.col(i)=Xai;
        }
    Xa		= mean(XaEn,1);
    Pa		= cov(trans(XaEn));
}

void EnsembleKalmanFilter::update(mat& XfEn)
{
    state_num  = XfEn.n_rows;
    ensemble_num  = XfEn.n_cols;
    Xa.set_size(state_num,1);
    Pa.set_size(state_num,state_num);
    XaEn.set_size(state_num,ensemble_num);
    Xa		= mean(XfEn,1);
    Pa		= cov(trans(XfEn));
    XaEn    = XfEn;
}

mat EnsembleKalmanFilter::ensembleInc(const mat& ensemble)
{
    unsigned int N	= ensemble.n_cols;
    mat	OneN;
    OneN.set_size(N, N);
    mat	I;			// Identity matrix (N*N)
    I.eye(N,N);

    OneN.fill(1.0/N); // 1N matrix (N*N)

    return ( ensemble*(I-OneN) );
}

void EnsembleKalmanFilter::perturbYo(const mat& m_R)
{
    mat ei;		// the ith member of observation perturbation matrix
    ei.set_size(observe_num, 1);

    mat mean1=zeros<mat>(observe_num, 1);
    Random yr;
    for (int i=0; i<ensemble_num; i++)
    {
        ei	= yr.MvNormal(mean1, m_R);
        YoEn.col(i)=Yo+ei;
    }
    R	= cov(trans(YoEn));	// R from samples
}

void EnsembleKalmanFilter::predict(const mat& m_Q)
{
    mat	Xfi(state_num, 1);		// the ith member of XfEn
    //Matrix	Q(m_Q);			// model error matrix
    mat	eta(state_num, 1);		// the ith member of model error
    mat	mean1=zeros<mat>(state_num, 1);
    //XfEn=ensemble;

    Random xr;

    for (int i=0; i<ensemble_num; i++)
    {
        eta	= xr.MvNormal(mean1, m_Q);
        Xfi = XfEn.col(i) + eta;

        //确保取值不超界
        if (Xl.n_rows==state_num && Xh.n_rows==state_num)
            for (int j=0; j<state_num; j++)
            {
                if ( Xfi(j, 0) < Xl(j, 0) )
                    Xfi(j,0)=Xl(j,0);
                if ( Xfi(j, 0) > Xh(j, 0) )
                    Xfi(j, 0)= Xh(j, 0);
            }

        XfEn.col(i)=Xfi;
    }

    XaEn = XfEn;	// where there is no update

    Xf		= mean(XfEn,1);
    Pf		= cov(trans(XfEn));		// not need! for clarification
}

mat EnsembleKalmanFilter::GetXa() const
{
    return Xa;
}

mat EnsembleKalmanFilter::GetXaEn() const
{
    return XaEn;
}

void EnsembleKalmanFilter::is_perturb_observe(const bool status)
{
    m_perturb_observe=status;
}
bool EnsembleKalmanFilter::is_perturb_observe()const
{
    return m_perturb_observe;
}

void EnsembleKalmanFilter::is_perturb_model(const bool status)
{
    m_perturb_model=status;
}
