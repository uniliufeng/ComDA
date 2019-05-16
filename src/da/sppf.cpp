#include "sppf.h"
using namespace ldas;

SPPF::SPPF()
{
    //ctor
}

SPPF::~SPPF()
{
    //dtor
}

SPPF::SPPF(const SPPF& other)
{
    //copy ctor
}

SPPF& SPPF::operator=(const SPPF& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

SPPF::SPPF(const itpp::mat particles_init, itpp::mat particlesNew_init, itpp::mat particlesPred_init, itpp::vec prior_init, itpp::Array < itpp::mat > SxPred_init, const itpp::vec obs,
           const float resampleThreshold)
{
    Xdim = particles_init.rows(); // extract state dimension
    Odim = obs.length(); // extract observation dimension
    U2dim = 1; // extract exogenous input 2 dimension
    Vdim = 1; // extract process noise dimension
    Ndim = 1; // extract observation noise dimension
    N = particles_init.cols(); // number of particles
    particles = particles_init; // copy particle buffer
    particlesNew = particlesNew_init;
    particlesPred = particlesPred_init;
    prior = prior_init;
    SxPred = SxPred_init;
    Obs = obs;
    weights = itpp::ones(N) / N; // particle weights
    St = itpp::round_i(resampleThreshold * N); // resample threshold
}

itpp::vec SPPF::GetXa() const
{
    return sum(elem_mult(reshape(repeat(weights, particles.rows()),
                                 particles.rows(), particles.cols()), particles), 2);;
}
itpp::mat SPPF::GetXaEn() const
{
    return particles;
}

///更新
void SPPF::update()
{
    itpp::vec normWeights(N);
    normWeights = 1.0 / N;

    itpp::mat UU2;
    UU2 = itpp::zeros(0, N);
    // EVALUATE IMPORTANCE WEIGHTS

    itpp::mat OBS = repeat(reshape(Obs, Odim, 1), N);

    // EVALUATE IMPORTANCE WEIGHTS

    // calculate observation likelihood for each particle (in log domain)
    itpp::vec likelihood_cal(N);
    likelihood_cal = likelihood(OBS, particlesPred) + 1e-99;

    itpp::mat difX;
    difX = particlesPred - particlesNew;

    itpp::vec proposal = itpp::zeros(N);
    double normfact = pow(2*PI_NUMBER,Xdim/2.0);

    for (int k = 0; k < N; k++)
    {
        itpp::mat cholFact = SxPred(k);
        itpp::vec foo = inv(cholFact) * difX.get_col(k);
        proposal(k) = exp(-0.5 * foo * foo) / abs(normfact * prod(
                          diag(cholFact))) + 1e-99;
        weights(k) = weights(k) * likelihood_cal(k) * prior(k) / proposal(k);
    }

    weights = weights / sum(weights);

    // RESAMPLE
    double S = 1.0 / sum(pow(weights, 2)); // calculate effective particle set size

    if (S < St)
    {
        itpp::ivec array(N);
        for (int j = 1; j <= N; j++)
        {
            array(j - 1) = j;
        }
        itpp::ivec outIndex = residualResample(array, weights);

        particles = particlesPred.get_cols(outIndex - 1);

//		for (int k = 0; k < N; k++) {
//			Sx(k) = SxPred(outIndex(k) - 1);
//		}
        weights = normWeights;
    }
    else
    {
        particles = particlesPred;
//		Sx = SxPred;
    }
}
