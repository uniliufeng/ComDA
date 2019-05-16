#include "srcdkf.h"
using namespace ldas;

SRCDKF::SRCDKF()
{
    //ctor
}

SRCDKF::~SRCDKF()
{
    //dtor
}

SRCDKF::SRCDKF(const SRCDKF& other)
{
    //copy ctor
}

SRCDKF& SRCDKF::operator=(const SRCDKF& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

SRCDKF::SRCDKF(const itpp::mat particles_init, const itpp::vec obs, const int xdim,
               const int odim, itpp::mat Z_, itpp::mat Sstate_Cov, itpp::mat Model_Cov, itpp::mat Obs_Cov)
{
    Xdim = xdim; // extract state dimension
    Odim = odim; // extract observation dimension
    U2dim = 1; // extract exogenous input 2 dimension
    Vdim = 1; // extract process noise dimension
    Ndim = 1; // extract observation noise dimension
    particles = particles_init; // copy particle buffer
    Obs = obs;
    cdkfParams = sqrt(3.0);
    Model_Noise_Mu = itpp::zeros(Xdim);
    Obs_Noise_Mu = itpp::zeros(Odim);
    Z = Z_;
    Sstate = Sstate_Cov; //lower triangular Cholesky factor of state covariance at time k-1
    Model_Covariance = Model_Cov;
    Obs_Covariance = Obs_Cov;
}

itpp::vec SRCDKF::GetXa() const
{
    return Xa;
}

///更新
void SRCDKF::update()
{

    // setup buffer
    itpp::vec xh = itpp::zeros(Xdim);
    itpp::vec xh_ = itpp::zeros(Xdim);
    itpp::vec yh_ = itpp::zeros(Odim);
    itpp::vec inov = itpp::zeros(Odim);

    // Get and calculate CDKF scaling parameters and sigma point weights
    double h = cdkfParams;
    double hh = itpp::sqr(h);

    int nsp1 = 2 * (Xdim + Vdim) + 1; // number of sigma points (first set)
    int nsp2 = 2 * (Xdim + Ndim) + 1; // number of sigma points (second set)

    // sigma-point weights set 1
    itpp::mat W1(2, 2); // sigma-point weights set 1
    W1(0, 0) = (hh - Xdim - Vdim) / hh;
    W1(0, 1) = 1 / (2 * hh);
    W1(1, 0) = 1 / (2 * h);
    W1(1, 1) = sqrt(hh - 1) / (2 * hh);

    itpp::mat W2 = W1; // sigma-point weights set 2
    W2(0, 0) = (hh - Xdim - Ndim) / hh;

    itpp::mat Zeros_Xdim_X_Vdim = itpp::zeros(Xdim, Vdim);
    itpp::mat Zeros_Vdim_X_Xdim = itpp::zeros(Vdim, Xdim);
    itpp::mat Zeros_Xdim_X_Ndim = itpp::zeros(Xdim, Ndim);
    itpp::mat Zeros_Ndim_X_Xdim = itpp::zeros(Ndim, Xdim);

    itpp::mat Sx, Sv, Sn;
    itpp::vec dv, dV, ds;
    itpp::ivec idx;
    double nu;
    int ind1 = 0, ind2 = 0, paramdim = 0;

    itpp::mat Sz(Xdim + Vdim, Xdim + Vdim);
    itpp::mat hSz(Xdim + Vdim, Xdim + Vdim);
    itpp::mat hSzM(Xdim + Vdim, Xdim + Vdim + Xdim + Vdim);

    itpp::mat Z2(Xdim + Ndim, nsp2);
    itpp::mat Sz2(Xdim + Ndim, Xdim + Ndim);
    itpp::mat hSz2(Xdim + Ndim, Xdim + Ndim);
    itpp::mat hSzM2(Xdim + Ndim, Xdim + Ndim + Xdim + Ndim);

    itpp::mat sub1Z, sub2Z, sub1Z2, sub2Z2;

    itpp::mat X_(Xdim, nsp1), A(Xdim, Xdim + Vdim), B(Xdim, Xdim + Vdim);
    itpp::mat C(Odim, Odim + Ndim), D(Odim, Odim + Ndim), Y_(Odim, nsp2);

    itpp::mat temp1, Sx_, AB, Sx1;
    itpp::mat temp2, Sy, Pxy(Xdim, Odim), KG(Xdim, Odim), CD;

    itpp::mat Syx1, Syw1, cov_update_vectors;
    itpp::mat temp3, compo;

    itpp::mat UU1, UU2;

    UU1 = itpp::zeros(0, nsp1);
    UU2 = itpp::zeros(0, nsp2);

    // TIME UPDATE
    Sx = Sstate;
    Sv = Model_Covariance; // get process noise covariance Cholesky factor
    Sn = Obs_Covariance; // get observation noise covariance Cholesky factor

    X_ = particles;

    xh_ = W1(0, 0) * X_.get_col(0) + W1(0, 1) * sum(
              X_(0, Xdim - 1, 1, nsp1 - 1), 2);
    A = W1(1, 0) * (X_(0, Xdim - 1, 1, Xdim + Vdim) - X_(0, Xdim - 1, Xdim
                    + Vdim + 1, nsp1 - 1));
    B = W1(1, 1) * (X_(0, Xdim - 1, 1, Xdim + Vdim) + X_(0, Xdim - 1, Xdim
                    + Vdim + 1, nsp1 - 1) - repeat(reshape(2 * X_.get_col(0), Xdim, 1),
                            Xdim + Vdim));

    AB = transpose(concat_horizontal(A, B));
    qr_economy(AB, temp1, Sx_);

    Sx_ = transpose(Sx_);
    Sx1 = Sx_;

    //--Calculate predicted observation mean, dealing with angular discontinuities if needed
    Z2 = repeat(reshape(concat(xh_, Obs_Noise_Mu), Xdim + Ndim, 1),
                nsp2);
    Sz2 = concat_vertical(concat_horizontal(Sx_, Zeros_Xdim_X_Ndim),
                          concat_horizontal(Zeros_Ndim_X_Xdim, Sn));
    hSz2 = h * Sz2;
    hSzM2 = concat_horizontal(hSz2, -hSz2);
    Z2.set_submatrix(0, Xdim + Ndim - 1, 1, nsp2 - 1, Z2(0, Xdim + Ndim - 1, 1,
                     nsp2 - 1) + hSzM2); // build sigma-point set

    sub1Z2 = Z2.get_rows(0, Xdim - 1);
    sub2Z2 = Z2.get_rows(Xdim, Xdim + Ndim - 1);

    Y_ = sub1Z2 + sub2Z2; // propagate through observation model

    yh_ = W2(0, 0) * Y_.get_col(0) + W2(0, 1) * sum(
              Y_(0, Odim - 1, 1, nsp2 - 1), 2);
    C = W2(1, 0) * (Y_(0, Odim - 1, 1, Xdim + Ndim) - Y_(0, Odim - 1, Xdim
                    + Ndim + 1, nsp2 - 1));
    D = W2(1, 1) * (Y_(0, Odim - 1, 1, Xdim + Ndim) + Y_(0, Odim - 1, Xdim
                    + Ndim + 1, nsp2 - 1) - repeat(reshape(2 * Y_.get_col(0), Odim, 1),
                            Xdim + Ndim));
    //-- Calculate predicted observation mean

    CD = transpose(concat_horizontal(C, D));

    qr_economy(CD, temp2, Sy);

    Sy = transpose(Sy);

    // MEASUREMENT UPDATE

    Syx1 = C(0, Odim - 1, 0, Xdim - 1);
    Syw1 = C(0, Odim - 1, Xdim, Xdim + Ndim - 1);

    Pxy = Sx_ * transpose(Syx1);

    KG = Pxy * inv(transpose(Sy)) * inv(Sy);

    inov = Obs - yh_; // inovation (observation error)
    xh = xh_ + KG * inov;

    compo = transpose(concat_horizontal(concat_horizontal(Sx_ - KG * Syx1, KG
                                        * Syw1), KG * D));
    qr_economy(compo, temp3, Sx);
    Sx = transpose(Sx);

    Xa = xh;
    //	Output.Sx = Sx;
    //	Output.Sx_ = Sx1;
    //	Output.InternalVariablesDS.xh_ = xh_;
    //	Output.InternalVariablesDS.Sx_ = Sx_;
    //	Output.InternalVariablesDS.yh_ = yh_;
    //	Output.InternalVariablesDS.inov = inov;
    //	Output.InternalVariablesDS.inov_cov = Sy;
    //	Output.InternalVariablesDS.KG = KG;
    //	Output.pNoise = InputDS.pNoise;
    //	Output.oNoise = InputDS.oNoise;

}
