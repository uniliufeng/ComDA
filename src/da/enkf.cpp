#include "enkf.hpp"
using namespace ldas;

EnKF::EnKF():m_seed(0),n(0),m(0),N(0),m_perturb_observe(true),m_perturb_model(true)
{

}

EnKF::~EnKF() {}

EnKF::EnKF(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum):m_seed(0),n(StateNum),N(EnNum),m(ObsNum),m_perturb_observe(true),m_perturb_model(true)
{
    Xf.setDim(n, 1);
    Xa.setDim(n, 1);
    XfEn.setDim(n, N);
    XaEn.setDim(n, N);
    Y.setDim(m, 1);
    YEn.setDim(m, N);
    Yo.setDim(m, 1);
    YoEn.setDim(m, N);
    R.setDim(m, m);
    Pf.setDim(n, n);
    Pa.setDim(n, n);
    K.setDim(n, m);
}

void EnKF::stateNum(const unsigned int num)
{
    n=num;
    Xf.setDim(n, 1);
    Xa.setDim(n, 1);
    Pf.setDim(n, n);
    Pa.setDim(n, n);
    if (m>0)
    {
        K.setDim(n, m);
    }
    if (N>0)
    {
        XfEn.setDim(n, N);
        XaEn.setDim(n, N);
    }
}
int EnKF::stateNum() const
{
    return n;
}

void EnKF::observeNum(const unsigned int num)
{
    m=num;
    Y.setDim(m, 1);
    Yo.setDim(m, 1);
    R.setDim(m, m);
    if (n>0)
    {
        K.setDim(n, m);
    }
    if (N>0)
    {
        YEn.setDim(m, N);
        YoEn.setDim(m, N);
    }
}
int EnKF::observeNum() const
{
    return m;
}

void EnKF::ensembleNum(const unsigned int num)
{
    N=num;
    if (n>0)
    {
        XfEn.setDim(n, N);
        XaEn.setDim(n, N);
    }
    if (m>0)
    {
        YEn.setDim(m, N);
        YoEn.setDim(m, N);
    }
}
int EnKF::ensembleNum() const
{
    return N;
}

unsigned long EnKF::randomSeed() const
{
    return m_seed;
}
void EnKF::randomSeed(const unsigned long seed)
{
    m_seed = seed;
}

void EnKF::is_perturb_observe(const bool status)
{
    m_perturb_observe=status;
}
bool EnKF::is_perturb_observe()const
{
    return m_perturb_observe;
}

void EnKF::is_perturb_model(const bool status)
{
    m_perturb_model=status;
}
bool EnKF::is_perturb_model()const
{
    return m_perturb_model;
}

void EnKF::xrange(const Matrix<double>& m_Xl, const Matrix<double>& m_Xh)
{
    Xl = m_Xl;
    Xh = m_Xh;
}

void EnKF::yrange(const Matrix<double>& m_Yl, const Matrix<double>& m_Yh)
{
    Yl = m_Yl;
    Yh = m_Yh;
}

void EnKF::observe(const Matrix<double>& m_Yo)
{
    this->Yo = m_Yo;
}

Matrix<double> EnKF::modelError(const Matrix<double>& m_Xl,const Matrix<double>& m_Xh, const double fac)
{
    Matrix<double> eta(n,1);		// vector of model error variance
    Matrix<double> Q(n, n);			// model error matrix

    eta		= fac*(m_Xh-m_Xl);
    Q		= eta*eta.tran();

    // make Q a positive definite matrix
    for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
            if ( i!=j )
                Q.set(i, j, 0.5*Q.get(i, j));

    return Q;
}

void EnKF::predict(const Matrix<double>& m_Q)
{
    Matrix<double>	Xfi(n, 1);		// the ith member of XfEn
    //Matrix	Q(m_Q);			// model error matrix
    Matrix<double>	eta(n, 1);		// the ith member of model error
    Matrix<double>	mean(n, 1);
    mean	= 0.0;
    //XfEn=ensemble;

    Random xr(m_seed);
    if (m_seed==0)
        xr=Random((unsigned)time(NULL));

    for (int i=1; i<=N; i++)
    {
        eta	= xr.MvNormal(mean, m_Q);
        Xfi = XfEn.extractColV(i) + eta;

        //确保取值不超界
        if (Xl.row()==n && Xh.row()==n)
            for (int j=1; j<=n; j++)
            {
                if ( Xfi.get(j, 1) < Xl.get(j, 1) )
                    Xfi.set(j, 1, Xl.get(j,1));
                if ( Xfi.get(j, 1) > Xh.get(j, 1) )
                    Xfi.set(j, 1, Xh.get(j, 1));
            }

        XfEn.setColV(i, Xfi);
    }

    XaEn = XfEn;	// where there is no update

    Xf		= ensembleMean(XfEn);
    Pf		= ensembleCov(XfEn);		// not need! for clarification
}

void EnKF::update(Matrix<double>& m_XfEn, Matrix<double>& m_YEn,
                  Matrix<double>& m_Yo, Matrix<double>& m_R, Matrix<double>& m_Q)
{
    n  = m_XfEn.row();
    N  = m_XfEn.column();
    m  = m_Yo.row();

    //Q=m_Q;
    XfEn = m_XfEn;
    YEn  = m_YEn;
    Yo   = m_Yo;
    R    = m_R;

    YoEn.setDim(m,N);
    Xa.setDim(n,1);
    Xf.setDim(n,1);
    XaEn.setDim(n,N);
    Pa.setDim(n,n);
    Pf.setDim(n,n);

    Matrix<double>	PH(n, m);			// Pf*H()
    Matrix<double>  HPH(m, m);          // H()*Pf*H()^T
    Matrix<double>	XfEnInc(n, N);		// Ensemble perturbation matrix of model forecast
    Matrix<double>	YEnInc(m, N);
    Matrix<double>	D(m, N);			// ensemble of inovation
    //Matrix<double>  D1(m,N);
    //Matrix<double>  variance(1,1);
    //Matrix<double>  Q1;
    Matrix<double>  temp(m,n);

    //generate XfEn with perturb model error
    if (m_perturb_model)
        predict(m_Q);

    //perturb observation
    if (m_perturb_observe)
        perturbYo(R);
    else
        for(int i=1; i<=N; i++)
            YoEn.setColV(i,Yo);
    D= YoEn - YEn;

    // increment
    XfEnInc	= ensembleInc(XfEn);
    YEnInc	= ensembleInc(YEn);
    Matrix<double> tmp(N,m);
    tmp=YEnInc.tran();

    PH		= XfEnInc*tmp/(N-1.0);
    HPH		= YEnInc *tmp/(N-1.0);
    temp = HPH+R;

    /*
    if(abs(temp.det())>1.0e-10)
    {
        K		= PH*temp.inv();
        XaEn	= XfEn + K*D;
    }
    else
    {
        XaEn    = XfEn;
    }*/
    try
    {
        K		= PH*temp.inv();
        XaEn	= XfEn + K*D;
    }
    catch (ldas::Exception e)
    {
        XaEn	= XfEn;
    }

    if (Xl.row()==n && Xh.row()==n)
        for (int i=1; i<=N; i++)
        {
            Matrix<double> Xai	= XaEn.extractColV(i);
            for (int j=1; j<=n; j++)
            {
                if ( Xai.get(j, 1) < Xl.get(j, 1) )
                    Xai.set(j, 1, Xl.get(j,1));
                if ( Xai.get(j, 1) > Xh.get(j, 1) )
                    Xai.set(j, 1, Xh.get(j, 1));
            }
            XaEn.setColV(i, Xai);
        }
    Xa		= ensembleMean(XaEn);
    Pa		= ensembleCov(XaEn);
}

void EnKF::update(Matrix<double>& XfEn)
{
    n  = XfEn.row();
    N  = XfEn.column();
    Xa.setDim(n,1);
    Pa.setDim(n,n);
    XaEn.setDim(n,N);
    Xa		= ensembleMean(XfEn);
    Pa		= ensembleCov(XfEn);
    XaEn    = XfEn;
}

Matrix<double> EnKF::GetXf() const
{
    return Xf;
}
Matrix<double> EnKF::GetXa() const
{
    return Xa;
}
Matrix<double> EnKF::GetXfEn() const
{
    return XfEn;
}
Matrix<double> EnKF::GetXaEn() const
{
    return XaEn;
}
Matrix<double> EnKF::GetY() const
{
    return Y;
}
Matrix<double> EnKF::GetYEn() const
{
    return YEn;
}
Matrix<double> EnKF::GetYo() const
{
    return Yo;
}
Matrix<double> EnKF::GetYoEn() const
{
    return YoEn;
}
Matrix<double> EnKF::GetR() const
{
    return R;
}
Matrix<double> EnKF::GetPf() const
{
    return Pf;
}
Matrix<double> EnKF::GetPa() const
{
    return Pa;
}
Matrix<double> EnKF::gain() const
{
    return K;
}

Matrix<double> EnKF::observeError(const Matrix<double>& m_Yl, const Matrix<double>& m_Yh, const double fac)
{
    Matrix<double> epsilon(m, 1);	// variance of observation error
    Matrix<double> noise(m, m);

    epsilon	= fac*(m_Yh-m_Yl);
    noise	= epsilon*epsilon.tran();
    // make noise a positive definite matrix
    for (int i=1; i<=m; i++)
    {
        for (int j=1; j<=m; j++)
        {
            if ( i!=j )
            {
                noise.set(i, j, 0.5*noise.get(i, j));
            }
        }
    }

    return noise;
}

void EnKF::perturbYo(const Matrix<double>& m_R)
{
    Matrix<double> ei(m, 1);		// the ith member of observation perturbation matrix

    Matrix<double> mean(m, 1);
    mean	= 0.0;
    Random yr(m_seed);
    if (m_seed==0)
        yr=Random((unsigned)time(NULL));
    for (int i=1; i<=N; i++)
    {
        ei	= yr.MvNormal(mean, m_R);
        // ei	= yr.TruncatedMvNormal(mean, noise, Yl, Yh);
        YoEn.setColV(i, Yo+ei);
    }
    R	= ensembleCov(YoEn);	// R from samples
}

Matrix<double> EnKF::ensembleInc(const Matrix<double>& ensemble)
{
    unsigned int N	= ensemble.column();
    Matrix<double>	OneN(N, N);			// 1N matrix (N*N)
    Matrix<double>	I(N, N);			// Identity matrix (N*N)

    OneN	= 1.0/(double)N;
    I.genI();

    return ( ensemble*(I-OneN) );
}

Matrix<double> EnKF::ensembleMean(const Matrix<double>& ensemble)
{
    unsigned int N	= ensemble.column();
    Matrix<double> mean(ensemble.row(), 1);

    for (int i=1; i<=N; i++)
    {
        mean += ensemble.extractColV(i);
    }
    mean /= N;

    return mean;
}

Matrix<double> EnKF::ensembleCov(const Matrix<double>& ensemble)
{
    Matrix<double> cov(ensemble.row(), ensemble.row());

    unsigned int N	= ensemble.column();
    Matrix<double>	OneN(N, N);			// 1N matrix (N*N)
    Matrix<double>	I(N, N);			// Identity matrix (N*N)

    OneN	= 1.0/(double)N;
    I.genI();

    cov = ensemble*(I-OneN)*ensemble.tran();
    cov /= (N-1.0);

    return cov;
}
