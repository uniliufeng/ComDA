#include "ekf.h"
using namespace ldas;

EKF::EKF()
{
    //ctor
}

EKF::~EKF()
{
    //dtor
}

EKF::EKF(const EKF& other)
{
    //copy ctor
}

EKF& EKF::operator=(const EKF& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

EKF::EKF(const itpp::vec xf_init, const itpp::vec obs, const int xdim, const int odim,
         itpp::mat C_, itpp::mat Sstate_Cov, itpp::mat Model_Cov, itpp::mat Obs_Cov)
{
    Xdim = xdim; // extract state dimension
    Odim = odim; // extract observation dimension
    U2dim = 1; // extract exogenous input 2 dimension
    Vdim = 1; // extract process noise dimension
    Ndim = 1; // extract observation noise dimension
    xf = xf_init; // copy particle buffer
    Obs = obs;
    Model_Noise_Mu = itpp::zeros(Xdim);
    Obs_Noise_Mu = itpp::zeros(Odim);
    C = C_;
    Sstate = Sstate_Cov; //lower triangular Cholesky factor of state covariance at time k-1
    Model_Covariance = Model_Cov;
    Obs_Covariance = Obs_Cov;
}

itpp::vec EKF::GetXa() const
{
    return Xa;
}

///更新
void EKF::update()
{

    itpp::mat Py = C * Px_ * transpose(C) + Obs_Covariance;
    itpp::mat KG = Px_ * transpose(C) * inv(Py);
    itpp::vec inov = Obs - Hx; // inovation (observation error)
    itpp::vec xh = xf + KG * inov;
    itpp::vec Xa = xh;
}
