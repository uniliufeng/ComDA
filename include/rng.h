#ifndef __LDAS_RNG_H
#define __LDAS_RNG_H

#include "random.hpp"
#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/itsignal.h>
#include <iostream>

namespace ldas
{
struct sn_struct
{
    itpp::vec xi;
    itpp::mat Omega;
    itpp::vec alpha;
    itpp::vec omega_rebel; //vector of scale parameters
    itpp::vec mean; //numeric vector representing the mean value of the distribution
    itpp::mat variance; //variance matrix of the distribution
    itpp::mat Omega_cor; //correlation matrix associated to Omega
    itpp::mat Omega_conc; //concentration matrix associated to Omega
    itpp::mat Psi; //covariance matrix of the equivalent (lambda,Psi) parametrization
    itpp::vec lambda; //shape parameters of the marginal distributions
    itpp::vec delta;
    itpp::vec skewness; //numeric vector with marginal indices of skewness
};

class Rng : public Random
{
public:
    /** Default constructor */
    Rng();
    /** Default destructor */
    virtual ~Rng();

    ///Generates a random Gaussian (0,1) vector.
    itpp::vec randn(int size);
    ///Generates a random Gaussian (0,1) matrix.
    itpp::mat randn(int rows, int cols);

    ///Generates a random uniform (0,1) vector.s
    itpp::vec randu(int size);
    ///Generates a random uniform (0,1) matrix.
    itpp::mat randu(int rows, int cols);

///Multivariate Skew-Normal random number generator:xi-mean Omega-covariace matrix
    itpp::mat rmsn(int n, itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha);
///Multivariate Skew-Normal PDF
    itpp::vec dmsn(itpp::mat &x, itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha);

    itpp::vec uniform_vec(const double a, const double b, int size);

    itpp::mat laplace_mat(itpp::vec &mu, itpp::vec &b, int rows, int cols);

/// This function returns a matrix associated with multivariate normal distribution with mean mu and covariance matrix sigma, cases is the number of samples
    itpp::mat mnormrnd(itpp::vec &mu, itpp::mat &sigma, int cases);
/// This function returns a vector associated with truncated multivariate normal distribution with mean mu and covariance matrix sigma, n is the number of samples
    itpp::vec tnorm_rnd(itpp::vec &amu, itpp::mat &sigma, itpp::vec &a, itpp::vec &b, itpp::vec &la, itpp::vec &lb);
///LHSNORM Generate a latin hypercube sample with a normal distribution
    itpp::mat lhsnorm(itpp::vec &mu, itpp::mat &sigma, int n);
    itpp::mat lhsnorm_trunctated(itpp::vec &mu, itpp::mat &sigma, int n, itpp::vec &a, itpp::vec &b, itpp::vec &la,
                                 itpp::vec &lb);

    /***** Subfunctions used by several functions *****/
    itpp::vec normcdf(itpp::vec &x);
    itpp::vec rank(itpp::vec &x);
    itpp::vec norminv(itpp::vec &p, double mu, double sigma);
    itpp::vec erfcinv(itpp::vec &y);

//Computes mean vector, variance matrix and other relevant quantities of a given multivariate skew-normal distribution.
    sn_struct msn_quantities(itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha);

protected:
private:
};
}
#endif // __LDAS_RNG_H
