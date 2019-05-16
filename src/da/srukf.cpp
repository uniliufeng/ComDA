#include "srukf.h"
using namespace ldas;

SRUKF::SRUKF()
{
    //ctor
}

SRUKF::~SRUKF()
{
    //dtor
}

SRUKF::SRUKF(const SRUKF& other)
{
    //copy ctor
}

SRUKF& SRUKF::operator=(const SRUKF& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

SRUKF::SRUKF(const itpp::mat particles_init, const itpp::vec obs, const int xdim, const int odim, itpp::mat Z_, itpp::mat Sstate_Cov,
             itpp::mat Model_Cov, itpp::mat Obs_Cov)
{
    Xdim = xdim; // extract state dimension
    Odim = odim; // extract observation dimension
    U2dim = 1; // extract exogenous input 2 dimension
    Vdim = 1; // extract process noise dimension
    Ndim = 1; // extract observation noise dimension
    particles = particles_init; // copy particle buffer
    Obs = obs;
    alpha = 1.0;
    beta = 2.0;
    kappa = 0.0;
    Model_Noise_Mu = itpp::zeros(Xdim);
    Obs_Noise_Mu = itpp::zeros(Odim);
    Z = Z_;
    Sstate = Sstate_Cov; //lower triangular Cholesky factor of state covariance at time k-1
    Model_Covariance = Model_Cov;
    Obs_Covariance = Obs_Cov;
}

itpp::vec SRUKF::GetXa() const
{
    return Xa;
}

itpp::mat SRUKF::GetSx() const
{
    return Sx_Out;
}

///更新
void SRUKF::update()
{

    // setup buffer
    itpp::vec xh = itpp::zeros(Xdim);
    itpp::vec xh_ = itpp::zeros(Xdim);
    itpp::vec yh_ = itpp::zeros(Odim);
    itpp::vec inov = itpp::zeros(Odim);

    itpp::mat Zeros_Xdim_X_Vdim = itpp::zeros(Xdim, Vdim);
    itpp::mat Zeros_Vdim_X_Xdim = itpp::zeros(Vdim, Xdim);
    itpp::mat Zeros_Xdim_X_Ndim = itpp::zeros(Xdim, Ndim);
    itpp::mat Zeros_Ndim_X_Xdim = itpp::zeros(Ndim, Xdim);
    itpp::mat Zeros_XdimVdim_X_Ndim = itpp::zeros(Xdim + Vdim, Ndim);
    itpp::mat Zeros_Ndim_X_XdimVdim = itpp::zeros(Ndim, Xdim + Vdim);

    itpp::mat SzT, Sz, sS, sSz, sSzM, sSM;
    itpp::mat subZ1, subZ2, subZ3;
    itpp::mat Sx, Sv, Sn;
    itpp::vec dv, dV, ds;
    itpp::ivec idx;
    double nu;
    int ind1 = 0, ind2 = 0, paramdim = 0;
    itpp::mat X_, X_bps, temp1, Y_;
    itpp::mat foo1, Sx_, temp1_sub, Sx1;
    itpp::mat temp2, Py, Pxy, KG;
    itpp::mat foo2, Sy, temp2_sub;
    itpp::mat UU1, UU2;

    int L = Xdim + Vdim + Ndim; // augmented state dimension
    int nsp = 2 * L + 1; // number of sigma-points
    kappa = itpp::sqr(alpha) * (L + kappa) - L; // compound scaling parameter

    itpp::vec W(3);
    W(0) = kappa;
    W(1) = 0.5;
    W(2) = 0.0;
    W /= (L + kappa); // sigma-point weights
    W(2) = W(0) + (1.0 - itpp::sqr(alpha)) + beta;

    itpp::vec sqrtW(W);
    int possitive_W3 = (W(2) > 0); // is zero'th covariance weight possitive?
    sqrtW.set_subvector(0, 1, sqrt(W(0, 1))); // square root weights
    sqrtW(2) = sqrt(abs(W(2)));

    double Sqrt_L_plus_kappa = sqrt(L + kappa);

    UU1 = itpp::zeros(0, nsp);
    UU2 = itpp::zeros(0, nsp);

    Sx = Sstate;
    Sv = Model_Covariance; // get process noise covariance Cholesky factor
    Sn = Obs_Covariance; // get observation noise covariance Cholesky factor

//	Z = repeat(reshape(concat(Xf, Model_Noise_Mu, Obs_Noise_Mu), L, 1), nsp);
//	SzT = concat_vertical(concat_horizontal(Sx, Zeros_Xdim_X_Vdim),
//			concat_horizontal(Zeros_Vdim_X_Xdim, Sv));
//	Sz = concat_vertical(concat_horizontal(SzT, Zeros_XdimVdim_X_Ndim),
//			concat_horizontal(Zeros_Ndim_X_XdimVdim, Sn));
//	sSz = Sqrt_L_plus_kappa * Sz;
//	sSzM = concat_horizontal(sSz, -sSz);
//	Z.set_submatrix(0, L - 1, 1, nsp - 1, Z(0, L - 1, 1, nsp - 1) + sSzM); // build sigma-point set
//
//	//-- Calculate predicted state mean
//
//	subZ1 = Z.get_rows(0, Xdim - 1);
//	subZ2 = Z.get_rows(Xdim, Xdim + Vdim - 1);
//
//	// propagate sigma-points through process model
//	int dim = subZ1.rows();
//	int nov = subZ1.cols();
//	mat new_state = zeros(dim, nov);
//
//	vec state_col, V_col, U1_col;
//	for (int i = 0; i < nov; i++) {
//		state_col = subZ1.get_col(i);
//		V_col = subZ2.get_col(i);
//		U1_col = UU1.get_col(i);
//		state_col = InferenceDS.modelfunc.ffun(InferenceDS.model, state_col,
//				V_col, U1_col);
//		new_state.set_col(i, state_col);
//	}

    X_ = particles;

    X_bps = X_;
    xh_ = W(0) * X_.get_col(0) + W(1) * sum(X_.get_cols(1, nsp - 1), 2);
    temp1 = X_ - repeat(reshape(xh_, Xdim, 1), nsp);

    // QR update of state Cholesky factor. NOTE: here Sx_ is the UPPER Cholesky factor (Matlab excentricity)

    temp1_sub = sqrtW(1) * transpose(temp1.get_cols(1, nsp - 1));

    qr_economy(temp1_sub, foo1, Sx_);

    itpp::vec temp1_col = sqrtW(2) * temp1.get_col(0);

    if (possitive_W3) // deal with possible negative zero'th covariance weight
    {
        const char plus = '+';
        Sx_ = cholupdate(Sx_, temp1_col, plus);
    }
    else
    {
        const char minus = '-';
        Sx_ = cholupdate(Sx_, temp1_col, minus); // NOTE: here Sx_ is the UPPER Cholesky factor (Matlab excentricity)
    }

    Sx1 = transpose(Sx_);

    subZ3 = Z.get_rows(Xdim + Vdim, Xdim + Vdim + Ndim - 1);
    Y_ = X_bps + subZ3;

    //-- Calculate predicted observation mean

    yh_ = W(0) * Y_.get_col(0) + W(1) * sum(Y_.get_cols(1, nsp - 1), 2);
    temp2 = Y_ - repeat(reshape(yh_, Odim, 1), nsp);

    // QR update of observation error Cholesky factor. NOTE: here Sy is the UPPER Cholesky factor (Matlab excentricity)

    temp2_sub = sqrtW(1) * transpose(temp2.get_cols(1, nsp - 1));
    qr_economy(temp2_sub, foo2, Sy);

    itpp::vec temp2_col = sqrtW(2) * temp2.get_col(0);
    if (possitive_W3) // deal with possible negative zero'th covariance weight
    {
        const char plus = '+';
        Sy = cholupdate(Sy, temp2_col, plus);
    }
    else
    {
        const char minus = '-';
        Sy = cholupdate(Sy, temp2_col, minus); // NOTE: here Sy  is the UPPER Cholesky factor (Matlab excentricity)
    }

    Sy = transpose(Sy); // We need the lower triangular Cholesky factor

    // MEASUREMENT UPDATE

    Pxy = W(2) * (temp1.get_cols(0, 0) * transpose(temp2.get_cols(0, 0)))
          + W(1) * (temp1.get_cols(1, nsp - 1) * transpose(temp2.get_cols(1,
                    nsp - 1)));

    KG = Pxy * inv(transpose(Sy)) * inv(Sy);

    inov = Obs - yh_; // inovation (observation error)
    xh = xh_ + KG * inov;

    itpp::mat cov_update_vectors = KG * Sy; // Correct covariance. This is equivalent to :  Px = Px_ - KG*Py*KG';

    for (int j = 0; j < Odim; j++)
    {
        itpp::vec cov_update_vectors_col = cov_update_vectors.get_col(j);
        const char minus = '-';
        //cout<<Sx_<<" "<<cov_update_vectors_col<<endl;
        Sx_ = cholupdate(Sx_, cov_update_vectors_col, minus);
        //cout<<Sx_<<endl;
    }

    Sx = transpose(Sx_);

    Xa = xh;
    Sx_Out = Sx;
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
