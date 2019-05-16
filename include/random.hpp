#ifndef __LDAS_RANDOM_H
#define __LDAS_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif
#include "dSFMT.h"
#ifdef __cplusplus
}
#endif

#include "constant.h"
#include "matrix.hpp"
#include <armadillo>

namespace ldas
{
/// 随机数类，生成满足不同分布条件的随机数
class Random
{
public:
    /** \brief 默认构造函数
    *
    * 以默认参数/随机数种子构建
    */
    Random();
    ///指定随机数种子
    Random(const unsigned long seed);
    ///析构函数
    virtual ~Random();
    /** Copy constructor
     *  \param other Object to copy from
     */
    Random(const Random& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Random operator=(const Random& other);
    ///指定随机数种子
    void seed(const unsigned long s);
    ///获取随机数种子
    unsigned long seed() const;

    /** \brief Normal Distribution 满足高斯分布（正态分布）的随机变量
    *
    * This function returns a Gaussian random variate, with mean mu and standard deviation sigma.

    * The exponential distribution has the form:

    * \f$p(x) dx = {1 \over \sqrt{2\pi\sigma^2}} \exp (-x^2 / 2\sigma^2) dx\f$

    * for x in the range \f$(-\infty,+\infty)\f$
    * \param mu 期望/均值
    * \param sigma 标准差
    * \return 随机数
    */
    double gaussian(const double mu, const double sigma);

    /** \brief Gaussian tail Distribution 截尾正态分布
    *
    * This function provides random variates from the upper tail of a Gaussian distribution with standard deviation sigma. The values returned are larger than the lower limit a, which must be positive

    The probability distribution for Gaussian tail random variates is,

    \f$ p(x) dx = {1 \over N(a;\sigma)} \exp (- x^2/(2 \sigma^2)) dx \f$

    for x > a where \f$N(a;\sigma)\f$ is the normalization constant,

    \f$ N(a;\sigma) = (1/2) erfc(a / sqrt(2 sigma^2))\f$.
    * \param a
    * \param sigma
    * \return 随机数
    */
    double gaussian_tail(const double a, const double sigma);

    ///Generates a random uniform (0,1) number.
    double randu();

    /** \brief The Bivariate Gaussian probability distribution is

    This function generates a pair of correlated gaussian variates, with mean zero, correlation coefficient rho and standard deviations sigma_x and sigma_y in the x and y directions

      \f$ p(x,y) dxdy = (1/(2 pi sigma_x sigma_y sqrt(r)))
                        exp(-(x^2 + y^2 - 2 r x y)/(2c)) dxdy \f$

       for x,y in the range \f$ (-\infty,+\infty) \f$.

       The correlation coefficient rho should lie between 1 and -1
    */
    void bivariate_gaussian(double sigma_x, double sigma_y, double rho, double *x, double *y);

    /** \todo normal cdf? */
    double normal_01_cdf(double x);

    double normal_01_cdf_inv(double p);

    ///This function returns a random variate from the exponential distribution with mean mu
    double exponential(const double mu);
    ///This function returns a random variate from the Laplace distribution with width a
    double laplace(const double a);
    ///This function returns a random variate from the Cauchy distribution with scale parameter a
    double cauchy(const double a);
    ///This function returns a random variate from the Rayleigh distribution with scale parameter sigma
    double rayleigh(const double sigma);
    ///This function returns a random variate from the tail of the Rayleigh distribution with scale parameter sigma and a lower limit of a
    double rayleigh_tail (const double a, const double sigma);
    ///This function returns a random variate from the gamma distribution
    double gamma(const double a, const double b);
    ///This function returns a random variate from the uniform(flat) distribution from a to b
    double uniform(const double a, const double b);
    ///This function returns a random variate from the lognormal distribution
    double lognormal(const double zeta, const double sigma);
    ///This function returns a random variate from the chi-squared distribution with nu degrees of freedom
    double chisq(const double nu);
    ///This function returns a random variate from the F-distribution with degrees of freedom nu1 and nu2
    double fdist(const double nu1, const double nu2);
    ///This function returns a random variate from the t-distribution
    double tdist (const double nu);
    ///This function returns a random variate from the beta distribution
    double beta(const double a, const double b);
    ///This function returns a random variate from the logistic distribution
    double logistic(const double a);
    ///This function returns a random variate from the Pareto distribution of order a
    double pareto(double a, const double b);
    ///This function returns a random variate from the Weibull distribution
    double weibull (const double a, const double b);
    ///This function returns a random variate from the Type-1 Gumbel distribution
    double gumbel1 (const double a, const double b);
    ///This function returns a random variate from the Type-2 Gumbel distribution
    double gumbel2(const double a, const double b);
    ///This function returns an array of K random variates from a Dirichlet distribution of order K-1
    void dirichlet(const unsigned long K, const double alpha[], double theta[]);
    ///This function returns a random integer from the Poisson distribution with mean mu
    unsigned long poisson(double mu);
    ///This function returns either 0 or 1, the result of a Bernoulli trial with probability p
    unsigned long bernoulli(double p);
    ///This function returns a random integer from the binomial distribution, the number of successes in n independent trials with probability p
    unsigned long binomial(double p, unsigned long n);
    ///This function returns an array of K random variates from a multinomial distribution
    ///void multinomial(const unsigned int K, const unsigned long N, const double p[], unsigned long n[]);
    ///This function returns a random integer from the negative binomial distribution, the number of failures occurring before n successes in independent trials with probability p of success
    unsigned long negative_binomial(double p, double n);
    ///This function returns a random integer from the Pascal distribution. The Pascal distribution is simply a negative binomial distribution with an integer value of n
    /// unsigned long pascal(double p, unsigned long n);
    ///This function returns a random integer from the geometric distribution, the number of independent trials with probability p until the first success
    unsigned long geometric(const double p);
    ///This function returns a random integer from the hypergeometric distribution
    unsigned long hypergeometric(unsigned long n1, unsigned long n2, unsigned long t);
    ///This function returns a random integer from the logarithmic distribution
    unsigned long logarithmic(const double p);

    ///univariate Skew-Normal random number generator: location--mean scale-standard deviation
    double rsn(double location, double scale, double shape);
    ///univariate Skew-Normal PDF
    double dsn(double x, double location, double scale, double shape);
	///Generates a random Gaussian (0,1) variable.
    double randn();
	/// Computes the probability density p(x) at x for a Gaussian distribution with mean mu and standard deviation sigma
    double gaussian_pdf(const double x, const double mu=0, const double sigma=1);
    /** \todo d_huge/dpoly_value? */
    double d_huge ( void );
    double dpoly_value ( int n, double a[], double x );

    /// This function returns n points associated with the n dimensional truncated normal distribution with mean mu and covariance matrix sigma
    Matrix<double> tnorm_rnd(int n, Matrix<double> amu, Matrix<double> sigma, Matrix<double> a, Matrix<double> b,Matrix<double> la, Matrix<double> lb,Matrix<double> d, Matrix<int> kstep);
    /// This function returns a random number from a normal truncated to (left,right) interval with mean = mu, variance = sigma2
    double normt_rnd(double mu,double sigma2,double left,double right);
    /// This function returns a random number from a left-truncated normal distribution, with mean = mu, variance = sigma2
    double normlt_rnd(double mu,double sigma2,double left);
    /// This function returns a random number from a right-truncated normal distribution, with mean = mu, variance = sigma2
    double normrt_rnd(double mu,double sigma2,double right);

    /// multi-variable truncated normal distribution
    Matrix<double> MvNormal(Matrix<double>& mean, Matrix<double>& covariance, Matrix<double>& MinValue, Matrix<double> & MaxValue);
    /// multi-variable truncated normal distribution
    Matrix<double> MvNormal(Matrix<double>& mean, Matrix<double>& covariance, Matrix<double>& MinValue, Matrix<double> & MaxValue, int SampleNum);
    /// multi-variable truncated normal distribution
    Matrix<double> MvNormal(Matrix<double>& mean, Matrix<double>& covariance, int SampleNum);
    /** \brief 矩阵变量的正态分布.
    * generates a vector of multivariate normal distribution
    */
    Matrix<double> MvNormal(const Matrix<double>& mean, const Matrix<double>& covariance);
    arma::mat MvNormal(const arma::mat& mean, const arma::mat& covariance);
protected:
    /* generates a random number on [0,1]-real-interval */
    double genrand_real1(void);

    /* generates a random number on [0,1)-real-interval */
    double genrand_real2(void);

    /* generates a random number on (0,1)-real-interval */
    double genrand_real3(void);
private:
    // private member functions used by several functions
    double ugaussian();
    // generates a random number on [min,max]-interval
    unsigned long genrand_int(const unsigned long min, const unsigned long max);
    double gamma_int(const unsigned long a);
    double gamma_large(const double a);
    double gamma_frac(const double a);
    ///import dSFMT support
    dsfmt_t dsfmt;
    unsigned long m_seed;
};
}
#endif // __LDAS_RANDOM_H
