#include "random.hpp"
#include <cmath>

using namespace std;
namespace ldas
{
Random::Random():m_seed(0)
{
    dsfmt_init_gen_rand(&dsfmt,m_seed);
}

Random::Random(const unsigned long seed):m_seed(seed)
{
    dsfmt_init_gen_rand(&dsfmt,m_seed);
}

Random::~Random()
{
}

Random::Random(const Random& other)
{
    m_seed=other.seed();
    dsfmt_init_gen_rand(&dsfmt,m_seed);
}

Random Random::operator=(const Random& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    m_seed=other.seed();
    dsfmt_init_gen_rand(&dsfmt,m_seed);
    return *this;
}

void Random::seed(const unsigned long s)
{
    m_seed=s;
    dsfmt_init_gen_rand(&dsfmt,m_seed);
}

unsigned long Random::seed() const
{
    return m_seed;
}

/* generates a random number on [min,max]-interval */
unsigned long Random::genrand_int(const unsigned long min,const unsigned long max)
{
    unsigned long r;
    r = (unsigned long)((max - min + 1) * genrand_real2 ()) + min; // multiply interval with random and truncate

    if (r > max)
    {
        r = max;
    }
    if (max < min)
    {
        return 0x80000000;
    }

    return r;
}

/* generates a random number on [0,1]-real-interval */
double Random::genrand_real1(void)
{
    return dsfmt_genrand_open_close(&dsfmt);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double Random::genrand_real2(void)
{
    return dsfmt_genrand_close_open(&dsfmt);
}

/* generates a random number on (0,1)-real-interval */
double Random::genrand_real3(void)
{
    return dsfmt_genrand_open_open(&dsfmt);
}

/* These real versions are due to Isaku Wada, 2002/01/09 added */

/**** The Normal Distribution ****/
/* The exponential distribution has the form
p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) dx

for x in the range -\infty to +\infty*/
double Random::gaussian(const double mu, const double sigma)
{
    double x, y, r2;
    // make two normally distributed variates by Box-Muller transformation
    do
    {
        /* choose x,y in uniform square (-1,-1) to (+1,+1) */
        x = -1 + 2 * genrand_real2 ();
        y = -1 + 2 * genrand_real2 ();
        /* see if it is in the unit circle */
        r2 = x * x + y * y;
    }
    while (r2 > 1.0 || r2 == 0);
    /* Box-Muller transform */
    return (mu + sigma * y * std::sqrt (-2.0 * std::log (r2) / r2));
}

double Random::normal_01_cdf( double x )
{
    double a1 = 0.398942280444E+00;
    double a2 = 0.399903438504E+00;
    double a3 = 5.75885480458E+00;
    double a4 = 29.8213557808E+00;
    double a5 = 2.62433121679E+00;
    double a6 = 48.6959930692E+00;
    double a7 = 5.92885724438E+00;
    double b0 = 0.398942280385E+00;
    double b1 = 3.8052E-08;
    double b2 = 1.00000615302E+00;
    double b3 = 3.98064794E-04;
    double b4 = 1.98615381364E+00;
    double b5 = 0.151679116635E+00;
    double b6 = 5.29330324926E+00;
    double b7 = 4.8385912808E+00;
    double b8 = 15.1508972451E+00;
    double b9 = 0.742380924027E+00;
    double b10 = 30.789933034E+00;
    double b11 = 3.99019417011E+00;
    double cdf;
    double q;
    double y;
    //
    //  |X| <= 1.28.
    //
    if ( fabs ( x ) <= 1.28 )
    {
        y = 0.5 * x * x;

        q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
                                 + a6 / ( y + a7 ) ) ) );
        //
        //  1.28 < |X| <= 12.7
        //
    }
    else if ( fabs ( x ) <= 12.7 )
    {
        y = 0.5 * x * x;

        q = std::exp ( - y ) * b0 / ( fabs ( x ) - b1
                                      + b2 / ( fabs ( x ) + b3
                                               + b4 / ( fabs ( x ) - b5
                                                       + b6 / ( fabs ( x ) + b7
                                                               - b8 / ( fabs ( x ) + b9
                                                                       + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
        //
        //  12.7 < |X|
        //
    }
    else
    {
        q = 0.0;
    }
    //
    //  Take account of negative X.
    //
    if ( x < 0.0 )
    {
        cdf = q;
    }
    else
    {
        cdf = 1.0 - q;
    }

    return cdf;
}

/**** The Gaussian tail Distribution ****/
/* The probability distribution for Gaussian tail random variates is,

  p(x) dx = {1 \over N(a;\sigma)} \exp (- x^2/(2 \sigma^2)) dx

	for x > a where N(a;\sigma) is the normalization constant,

	  N(a;\sigma) = (1/2) erfc(a / sqrt(2 sigma^2)).
*/
double Random::gaussian_tail(const double a, const double sigma)
{
    /* Returns a gaussian random variable larger than a
    * This implementation does one-sided upper-tailed deviates.
    	*/

    double s = a / sigma;

    if (s < 1)
    {
        /* For small s, use a direct rejection method. The limit s < 1
        	can be adjusted to optimise the overall efficiency */

        double x;

        do
        {
            x = gaussian (0,1.0);
        }
        while (x < s);

        return x * sigma;
    }
    else
    {
        /* Use the "supertail" deviates from the last two steps
        * of Marsaglia's rectangle-wedge-tail method, as described
        * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
        * and the solution, p586.)
        	*/

        double u, v, x;

        do
        {
            u = genrand_real2 ();
            do
            {
                v = genrand_real2 ();
            }
            while (v == 0.0);
            x = std::sqrt (s * s - 2 * std::log (v));
        }
        while (x * u > s);

        return x * sigma;
    }
}

/**** The Bivariate Gaussian Distribution ****/
/* The Bivariate Gaussian probability distribution is

  p(x,y) dxdy = (1/(2 pi sigma_x sigma_y sqrt(r)))
  exp(-(x^2 + y^2 - 2 r x y)/(2c)) dxdy
  for x,y in the range -\infty to +\infty.
  The correlation coefficient rho should lie between 1 and -1

*/
void Random::bivariate_gaussian(double sigma_x, double sigma_y, double rho, double *x, double *y)
{
    double u, v, r2, scale;

    do
    {
        /* choose x,y in uniform square (-1,-1) to (+1,+1) */

        u = -1 + 2 * genrand_real2 ();
        v = -1 + 2 * genrand_real2 ();

        /* see if it is in the unit circle */
        r2 = u * u + v * v;
    }
    while (r2 > 1.0 || r2 == 0);

    scale = std::sqrt (-2.0 * std::log (r2) / r2);

    *x = sigma_x * u * scale;
    *y = sigma_y * (rho * u + std::sqrt(1 - rho*rho) * v) * scale;
}

Matrix<double> Random::MvNormal(const Matrix<double>& mean, const Matrix<double>& covariance)
{
    Matrix<double> vector(mean.row(), 1);
    Matrix<double> dec = covariance.chold();

    int i;
    for (i=1; i<=(int)mean.row(); i++)
    {
        vector.set(i, 1, gaussian(0,1));
    }
    vector = dec*vector + mean;

    return vector;
}

arma::mat Random::MvNormal(const arma::mat& mean, const arma::mat& covariance)
{
    arma::mat vector;
    vector.set_size(mean.n_rows, 1);
    arma::mat dec = chol(covariance);

    for (int i=0; i<mean.n_rows; i++)
    {
        vector(i, 0)=gaussian(0,1);
    }
    vector = dec*vector + mean;

    return vector;
}

/**** The Exponential Distribution ****/
/* The exponential distribution has the form

  p(x) dx = exp(-x/mu) dx/mu

for x = 0 ... +infty */
double Random::exponential(const double mu)
{
    double u = genrand_real3 ();
    return -mu * std::log (u);
}

/**** The Laplace Distribution ****/
/* The two-sided exponential probability distribution is

  p(x) dx = (1/(2 a)) * exp(-|x/a|) dx

for -infty < x < infty. It is also known as the Laplace distribution.  */
double Random::laplace(const double a)
{
    double u;
    do
    {
        u = 2 * genrand_real2 () - 1.0;
    }
    while (u == 0.0);

    if (u < 0)
    {
        return a * std::log (-u);
    }
    else
    {
        return -a * std::log (u);
    }
}

/**** The Cauchy Distribution ****/
/* The Cauchy probability distribution is

  p(x) dx = (1/(pi a)) (1 + (x/a)^2)^(-1) dx
  for x in the range -\infty to +\infty.
It is also known as the Lorentzian probability distribution */
double Random::cauchy (const double a)
{
    double u;
    do
    {
        u = genrand_real2 ();
    }
    while (u == 0.5);

    return a * std::tan (PI_NUMBER * u);
}

/**** The Rayleigh Distribution ****/
/* The Rayleigh distribution has the form

  p(x) dx = (x / sigma^2) exp(-x^2/(2 sigma^2)) dx

for x = 0 ... +infty */

double Random::rayleigh(const double sigma)
{
    double u = genrand_real3 ();

    return sigma * std::sqrt(-2.0 * std::log (u));
}

/**** The Rayleigh tailDistribution ****/
/* The Rayleigh tail distribution has the form

  p(x) dx = (x / sigma^2) exp((a^2 - x^2)/(2 sigma^2)) dx

for x = a ... +infty */

double Random::rayleigh_tail(const double a, const double sigma)
{
    double u = genrand_real3 ();

    return std::sqrt(a * a - 2.0 * sigma * sigma * std::log (u));
}

/**** The Gamma Distribution ****/
/* The Gamma distribution of order a>0 is defined by:

  p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

	for x>0.  If X and Y are independent gamma-distributed random
	variables of order a1 and a2 with the same scale parameter b, then
	X+Y has gamma distribution of order a1+a2.

The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */
double Random::gamma(const double a, const double b)
{
    /* assume a > 0 */
    unsigned long na = floor (a);

    if (a == na)
    {
        return b * genrand_int (0, na - 1);
    }
    else if (na == 0)
    {
        return b * gamma_frac (a);
    }
    else
    {
        return b * (genrand_int (0, na - 1) + gamma_frac (a - na)) ;
    }
}

/**** The Uniform Distribution ****/
/* This is the uniform distribution in the range [a, b)

  p(x) dx = 1/(b-a) dx   if  a <= x < b
  .....   = 0            otherwise

*/
double Random::uniform(const double a, const double b)
{
    double u = genrand_real2 ();

    /* A uniform distribution over [a,b] */

    return a * (1 - u) + b * u;
}

/**** The Lognormal Distribution ****/
/* The lognormal distribution has the form

  p(x) dx = 1/(x * sqrt(2 pi sigma^2)) exp(-(ln(x) - zeta)^2/2 sigma^2) dx

	for x > 0. Lognormal random numbers are the exponentials of
gaussian random numbers */
double Random::lognormal(const double zeta, const double sigma)
{
    double u, v, r2, normal, z;

    do
    {
        /* choose x,y in uniform square (-1,-1) to (+1,+1) */

        u = -1 + 2 * genrand_real2 ();
        v = -1 + 2 * genrand_real2 ();

        /* see if it is in the unit circle */
        r2 = u * u + v * v;
    }
    while (r2 > 1.0 || r2 == 0);

    normal = u * std::sqrt (-2.0 * std::log (r2) / r2);

    z =  std::exp (sigma * normal + zeta);

    return z;
}

/**** The Chi-squared Distribution ****/
/* The chi-squared distribution arises in statistics.
* If Y_i are n independent gaussian random variates with unit variance then the sum-of-squares,
* X_i = \sum_i Y_i^2
* has a chi-squared distribution with n degrees of freedom.
* The chisq distribution has the form
* p(x) dx = (1/(2*Gamma(nu/2))) (x/2)^(nu/2 - 1) exp(-x/2) dx
* for x = 0 ... +infty */
double Random::chisq(const double nu)
{
    double chisq = 2 * gamma (nu / 2, 1.0);
    return chisq;
}

/**** The F Distribution ****/
/* The F-distribution arises in statistics.
* If Y_1 and Y_2 are chi-squared deviates with \nu_1 and \nu_2 degrees of freedom then the ratio,
* X = { (Y_1 / \nu_1) \over (Y_2 / \nu_2) }
* has an F-distribution F(x;\nu_1,\nu_2).
* The F distribution has the form
* p(x) dx = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
Gamma(nu1/2) Gamma(nu2/2)) *
x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2) dx
The method used here is the one described in Knuth */

double Random::fdist(const double nu1, const double nu2)
{

    double Y1 =  gamma (nu1 / 2, 2.0);
    double Y2 =  gamma (nu2 / 2, 2.0);

    double f = (Y1 * nu2) / (Y2 * nu1);

    return f;
}

/**** The T Distribution ****/
/* The t-distribution arises in statistics.
* If Y_1 has a normal distribution and Y_2 has a chi-squared distribution with \nu degrees of freedom then the ratio,
* X = { Y_1 \over \sqrt{Y_2 / \nu} }
* has a t-distribution t(x;\nu) with \nu degrees of freedom.
* The t-distribution has the form

  p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
  * (1 + (x^2)/nu)^-((nu + 1)/2) dx
The method used here is the one described in Knuth */
double Random::tdist (const double nu)
{
    if (nu <= 2)
    {
        double Y1 = ugaussian ();
        double Y2 = chisq (nu);

        double t = Y1 / std::sqrt (Y2 / nu);

        return t;
    }
    else
    {
        double Y1, Y2, Z, t;
        do
        {
            Y1 = ugaussian ();
            Y2 = exponential (1 / (nu/2 - 1));

            Z = Y1 * Y1 / (nu - 2);
        }
        while (1 - Z < 0 || std::exp (-Y2 - Z) > (1 - Z));

        /* Note that there is a typo in Knuth's formula, the line below
        is taken from the original paper of Marsaglia, Mathematics of
        Computation, 34 (1980), p 234-256 */

        t = Y1 / std::sqrt ((1 - 2 / nu) * (1 - Z));
        return t;
    }
}

/**** The Beta Distribution ****/
/* The beta distribution has the form

  p(x) dx = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1) dx

The method used here is the one described in Knuth */
double Random::beta(const double a, const double b)
{
    double x1 = gamma (a, 1.0);
    double x2 = gamma (b, 1.0);

    return x1 / (x1 + x2);
}

/**** The Logistic Distribution ****/
/* The logistic distribution has the form,

  p(x) dx = (1/a) exp(-x/a) / (1 + exp(-x/a))^2 dx

for -infty < x < infty */
double Random::logistic(const double a)
{
    double x, z;

    do
    {
        x = genrand_real3 ();
    }
    while (x == 1);

    z = std::log (x / (1 - x));

    return a * z;
}

/**** The Pareto Distribution ****/
/* The Pareto distribution has the form,

  p(x) dx = (a/b) / (x/b)^(a+1) dx     for x >= b

*/
double Random::pareto(double a, const double b)
{
    double x = genrand_real3 ();

    double z = pow (x, -1 / a);

    return b * z;
}

/**** The Weibull Distribution ****/
/* The Weibull distribution has the form,

  p(x) dx = (b/a) (x/a)^(b-1) exp(-(x/a)^b) dx

*/
double Random::weibull(const double a, const double b)
{
    double x = genrand_real3 ();

    double z = pow (-std::log (x), 1 / b);

    return a * z;
}

/**** The Type I Gumbel Distribution ****/
/* The Type I Gumbel distribution has the form,

  p(x) dx = a b exp(-(b exp(-ax) + ax)) dx
  for -\infty < x < \infty.
*/
double Random::gumbel1(const double a, const double b)
{
    double x = genrand_real3 ();

    double z = (std::log(b) - std::log(-std::log(x))) / a;

    return z;
}

/**** The Type II Gumbel Distribution ****/
/* The Type II Gumbel distribution has the form,

  p(x) dx = b a x^-(a+1) exp(-b x^-a)) dx
  for 0 < x < \infty.
*/
double Random::gumbel2(const double a, const double b)
{
    double x = genrand_real3 ();

    double z = pow(-b / std::log(x), 1/a);

    return z;
}

/**** The Dirichlet probability Distribution ****/
/* The Dirichlet probability distribution of order K-1 is

  p(\theta_1,...,\theta_K) d\theta_1 ... d\theta_K =
  (1/Z) \prod_i=1,K \theta_i^{alpha_i - 1} \delta(1 -\sum_i=1,K \theta_i)

	The normalization factor Z can be expressed in terms of gamma functions:

      Z = {\prod_i=1,K \Gamma(\alpha_i)} / {\Gamma( \sum_i=1,K \alpha_i)}

		The K constants, \alpha_1,...,\alpha_K, must be positive. The K parameters,
		\theta_1,...,\theta_K are nonnegative and sum to 1.

		  The random variates are generated by sampling K values from gamma
		  distributions with parameters a=\alpha_i, b=1, and renormalizing.
		  See A.M. Law, W.D. Kelton, Simulation Modeling and Analysis (1991).

			Gavin E. Crooks <gec@compbio.berkeley.edu> (2002)
*/

void Random::dirichlet(const unsigned long K, const double alpha[], double theta[])
{
    unsigned long i;
    double norm = 0.0;

    for (i = 0; i < K; i++)
    {
        theta[i] = gamma (alpha[i], 1.0);
    }

    for (i = 0; i < K; i++)
    {
        norm += theta[i];
    }

    for (i = 0; i < K; i++)
    {
        theta[i] /= norm;
    }
}

/**** The Poisson Distribution ****/
/* The poisson distribution has the form

  p(n) = (mu^n / n!) exp(-mu)

for n = 0, 1, 2, ... . The method used here is the one from Knuth. */
unsigned long Random::poisson(double mu)
{
    double emu;
    double prod = 1.0;
    unsigned long k = 0;

    while (mu > 10)
    {
        unsigned long m = mu * (7.0 / 8.0);

        double X = gamma_int (m);

        if (X >= mu)
        {
            return k + binomial (mu / X, m - 1);
        }
        else
        {
            k += m;
            mu -= X;
        }
    }

    /* This following method works well when mu is small */

    emu = std::exp (-mu);

    do
    {
        prod *= genrand_real2 ();
        k++;
    }
    while (prod > emu);

    return k - 1;
}

/**** The Bernoulli Distribution ****/
/* The bernoulli distribution has the form,

  prob(0) = 1-p, prob(1) = p

*/
unsigned long Random::bernoulli(double p)
{
    double u = genrand_real2 () ;

    if (u < p)
    {
        return 1 ;
    }
    else
    {
        return 0 ;
    }
}

/**** The Binomial Distribution ****/
/* The binomial distribution has the form,

  prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n

This is the algorithm from Knuth */
unsigned long Random::binomial(double p, unsigned long n)
{
    unsigned long i, a, b, k = 0;

    while (n > 10)        /* This parameter is tunable */
    {
        double X;
        a = 1 + (n / 2);
        b = 1 + n - a;

        X = beta ((double) a, (double) b);

        if (X >= p)
        {
            n = a - 1;
            p /= X;
        }
        else
        {
            k += a;
            n = b - 1;
            p = (p - X) / (1 - X);
        }
    }

    for (i = 0; i < n; i++)
    {
        double u = genrand_real2 ();
        if (u < p)
        {
            k++;
        }
    }

    return k;
}

/**** The Multinomial Distribution ****/
/* The multinomial distribution has the form

  N!           n_1  n_2      n_K
  prob(n_1, n_2, ... n_K) = -------------------- p_1  p_2  ... p_K
  (n_1! n_2! ... n_K!)

	where n_1, n_2, ... n_K are nonnegative integers, sum_{k=1,K} n_k = N,
	and p = (p_1, p_2, ..., p_K) is a probability distribution.

	  Random variates are generated using the conditional binomial method.
	  This scales well with N and does not require a setup step.

		Ref:
		C.S. David, The computer generation of multinomial random variates,
		Comp. Stat. Data Anal. 16 (1993) 205-217
*/

//void Random::multinomial(const unsigned int K, const unsigned long N, const double p[], unsigned long n[])
//{
//	unsigned int k;
//	double norm = 0.0;
//	double sum_p = 0.0;
//
//	unsigned int sum_n = 0;
//
//  /* p[k] may contain non-negative weights that do not sum to 1.0.
//   * Even a probability distribution will not exactly sum to 1.0
//   * due to rounding errors.
//   */
//
//	for (k = 0; k < K; k++)
//    {
//		norm += p[k];
//    }
//
//	for (k = 0; k < K; k++)
//    {
//		if (p[k] > 0.0)
//        {
//			n[k] = binomial (p[k] / (norm - sum_p), N - sum_n);
//        }
//		else
//        {
//			n[k] = 0;
//        }
//
//		sum_p += p[k];
//		sum_n += n[k];
//    }
//
//}

/**** The Negative binomial Distribution ****/
/* The negative binomial distribution has the form,

  prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k

	for k = 0, 1, ... . Note that n does not have to be an integer.

This is the Leger's algorithm (given in the answers in Knuth) */
unsigned long Random::negative_binomial(double p, double n)
{
    double X = gamma (n, 1.0) ;
    unsigned long k = poisson (X*(1-p)/p) ;
    return k ;
}

/**** The Pascal Distribution ****/
/* The Pascal distribution is a negative binomial with valued integer n

  prob(k) =  (n - 1 + k)!/(n!(k - 1)!) *  p^n (1-p)^k for k = 0, 1, ..., n

*/
/*
unsigned long Random::pascal(double p, unsigned long n)
{
// This is a separate interface for the pascal distribution so that
// it can be optimized differently from the negative binomial in
// future.

  // e.g. if n < 10 it might be faster to generate the Pascal
  // distributions as the sum of geometric variates directly.

	unsigned long k = negative_binomial (p, (double) n);
	return k;
	}
*/

/**** The Geometric Distribution ****/
/* Geometric distribution (bernoulli trial with probability p)

  prob(k) =  p (1 - p)^(k-1) for n = 1, 2, 3, ...

	It gives the distribution of "waiting times" for an event that
occurs with probability p. */
unsigned long Random::geometric(const double p)
{
    double u = genrand_real3 ();

    unsigned long k;

    if (p == 1)
    {
        k = 1;
    }
    else
    {
        k = std::log (u) / std::log (1 - p) + 1;
    }

    return k;
}

/**** The Hypergeometric Distribution ****/
/* The hypergeometric distribution has the form,

  prob(k) =  choose(n1,t) choose(n2, t-k) / choose(n1+n2,t)

	where choose(a,b) = a!/(b!(a-b)!)

	  n1 + n2 is the total population (tagged plus untagged)
	  n1 is the tagged population
	  t is the number of samples taken (without replacement)
	  k is the number of tagged samples found

*/

unsigned long Random::hypergeometric(unsigned long n1, unsigned long n2, unsigned long t)
{
    const unsigned int n = n1 + n2;

    unsigned int i = 0;
    unsigned int a = n1;
    unsigned int b = n1 + n2;
    unsigned int k = 0;

    if (t > n)
    {
        t = n ;
    }

    if (t < n / 2)
    {
        for (i = 0 ; i < t ; i++)
        {
            double u = genrand_real2 () ;

            if (b * u < a)
            {
                k++ ;
                if (k == n1)
                    return k ;
                a-- ;
            }
            b-- ;
        }
        return k;
    }
    else
    {
        for (i = 0 ; i < n - t ; i++)
        {
            double u = genrand_real2 () ;

            if (b * u < a)
            {
                k++ ;
                if (k == n1)
                    return n1 - k ;
                a-- ;
            }
            b-- ;
        }
        return n1 - k;
    }
}

/**** The Logarithmic Distribution ****/
/* Logarithmic distribution

  prob(n) =   p^n / (n log(1/(1-p)) for n = 1, 2, 3, ...

	We use Kemp's second accelerated generator, from Luc Devroye's book
on "Non-Uniform Random Variate Generation", Springer */

unsigned long Random::logarithmic(const double p)
{
    double c = std::log (1-p) ;

    double v = genrand_real3 ();

    if (v >= p)
    {
        return 1 ;
    }
    else
    {
        double u = genrand_real3 ();
        double q = 1 - std::exp (c * u);

        if (v <= q*q)
        {
            double x = 1 + std::log(v)/std::log(q) ;
            return x ;
        }
        else if (v <= q)
        {
            return 2;
        }
        else
        {
            return 1 ;
        }
    }
}

/* Others  */
double Random::ugaussian()
{
    return gaussian (0,1.0);
}

double Random::gamma_int(const unsigned long a)
{
    if (a < 12)
    {
        unsigned long i;
        double prod = 1;

        for (i = 0; i < a; i++)
        {
            prod *= genrand_real3 ();
        }

        /* Note: for 12 iterations we are safe against underflow, since
        the smallest positive random number is O(2^-32). This means
        the smallest possible product is 2^(-12*32) = 10^-116 which
        is within the range of double precision. */

        return -std::log (prod);
    }
    else
    {
        return gamma_large((double) a);
    }
}

double Random::gamma_large(const double a)
{
    /* Works only if a > 1, and is most efficient if a is large

      This algorithm, reported in Knuth, is attributed to Ahrens.  A
      faster one, we are told, can be found in: J. H. Ahrens and
    	U. Dieter, Computing 12 (1974) 223-246.  */

    double sqa, x, y, v;
    sqa = std::sqrt (2 * a - 1);
    do
    {
        do
        {
            y = std::tan (PI_NUMBER * genrand_real2 ());
            x = sqa * y + a - 1;
        }
        while (x <= 0);
        v = genrand_real2 ();
    }
    while (v > (1 + y * y) * std::exp ((a - 1) * std::log (x / (a - 1)) - sqa * y));

    return x;
}

double Random::gamma_frac(const double a)
{
    /* This is exercise 16 from Knuth; see page 135, and the solution is
    	on page 551.  */

    double p, q, x, u, v;
    p = M_E / (a + M_E);
    do
    {
        u = genrand_real2 ();
        v = genrand_real3 ();

        if (u < p)
        {
            x = std::exp ((1 / a) * std::log (v));
            q = std::exp (-x);
        }
        else
        {
            x = 1 - std::log (v);
            q = std::exp ((a - 1) * std::log (x));
        }
    }
    while (genrand_real2 () >= q);

    return x;
}

double Random::d_huge ( void )
{
    return HUGE_VAL;
}

double Random::dpoly_value ( int n, double a[], double x )
{
    int i;
    double value;

    value = 0.0;

    for ( i = n-1; 0 <= i; i-- )
    {
        value = value * x + a[i];
    }

    return value;
}

double Random::normal_01_cdf_inv(double p)
{
    double a[8] = {3.3871328727963666080,     1.3314166789178437745e+2,
                   1.9715909503065514427e+3,  1.3731693765509461125e+4,
                   4.5921953931549871457e+4,  6.7265770927008700853e+4,
                   3.3430575583588128105e+4,  2.5090809287301226727e+3
                  };
    double b[8] = {	1.0,                       4.2313330701600911252e+1,
                    6.8718700749205790830e+2,  5.3941960214247511077e+3,
                    2.1213794301586595867e+4,  3.9307895800092710610e+4,
                    2.8729085735721942674e+4,  5.2264952788528545610e+3
                  };
    double c[8] = {	1.42343711074968357734,     4.63033784615654529590,
                    5.76949722146069140550,     3.64784832476320460504,
                    1.27045825245236838258,     2.41780725177450611770e-1,
                    2.27238449892691845833e-2,  7.74545014278341407640e-4
                  };
    double const1 = 0.180625;
    double const2 = 1.6;
    double d[8] = {	1.0,                        2.05319162663775882187,
                    1.67638483018380384940,     6.89767334985100004550e-1,
                    1.48103976427480074590e-1,  1.51986665636164571966e-2,
                    5.47593808499534494600e-4,  1.05075007164441684324e-9
                  };
    double e[8] = { 6.65790464350110377720,     5.46378491116411436990,
                    1.78482653991729133580,     2.96560571828504891230e-1,
                    2.65321895265761230930e-2,  1.24266094738807843860e-3,
                    2.71155556874348757815e-5,  2.01033439929228813265e-7
                  };
    double f[8] = { 1.0,                        5.99832206555887937690e-1,
                    1.36929880922735805310e-1,  1.48753612908506148525e-2,
                    7.86869131145613259100e-4,  1.84631831751005468180e-5,
                    1.42151175831644588870e-7,  2.04426310338993978564e-15
                  };
    double q;
    double r;
    double split1 = 0.425;
    double split2 = 5.0;
    double value;
    if ( p <= 0.0 )
    {
        value = -d_huge ( );
        return value;
    }
    if ( 1.0 <= p )
    {
        value = d_huge ( );
        return value;
    }
    q = p - 0.5;
    if ( fabs ( q ) <= split1 )
    {
        r = const1 - q * q;
        value = q * dpoly_value ( 8, a, r ) / dpoly_value ( 8, b, r );
    }
    else
    {
        if ( q < 0.0 )
        {
            r = p;
        }
        else
        {
            r = 1.0 - p;
        }
        if ( r <= 0.0 )
        {
            value = -1.0;
            return 1;
        }
        r = std::sqrt ( -std::log ( r ) );
        if ( r <= split2 )
        {
            r = r - const2;
            value = dpoly_value ( 8, c, r ) / dpoly_value ( 8, d, r );
        }
        else
        {
            r = r - split2;
            value = dpoly_value ( 8, e, r ) / dpoly_value ( 8, f, r );
        }
        if ( q < 0.0 )
        {
            value = -value;
        }
    }
    return value;
}

/**** The Truncated Gaussian Distribution ****/
double Random::normt_rnd(double mu,double sigma2,double left,double right)
{
    double std, lowerProb, upperProb, u, result;

    std = std::sqrt(sigma2);

    //Calculate bounds on probabilities
    lowerProb = normal_01_cdf((left-mu)/std);
    upperProb = normal_01_cdf((right-mu)/std);
    //Draw uniform from within (lowerProb,upperProb)
    u = lowerProb+(upperProb-lowerProb) * genrand_real3();
    //Find needed quantiles
    result = mu + normal_01_cdf_inv(u) * std;
    return result;
}
/**** The Left Truncated Gaussian Distribution ****/
double Random::normlt_rnd(double mu,double sigma2,double left)
{
    double right, std, result;
    std = std::sqrt(sigma2);
    right = mu - 5 * std;

    result = normt_rnd(mu,sigma2,left,right);

    if (result < left)
    {
        return left;
    }
    else
    {
        return result;
    }
}
/**** The Right Truncated Gaussian Distribution ****/
double Random::normrt_rnd(double mu,double sigma2, double right)
{
    double left, std, result;
    std = std::sqrt(sigma2);
    left = mu - 5 * std;

    return result = normt_rnd(mu,sigma2,left,right);
}

// This function returns n points associated with the n dimensional truncated normal distribution with mean mu and covariance matrix sigma
Matrix<double> Random::tnorm_rnd(int n,            Matrix<double> amu, Matrix<double> sigma,
                                 Matrix<double> a, Matrix<double> b,   Matrix<double> la,
                                 Matrix<double> lb,Matrix<double> d,   Matrix<int> kstep)
{
    int niter=10;
    double t1,aa;
    //transform to work in terms of z=d*x
    Matrix<double> z(n,1);
    Matrix<double> anu(n,1);
    Matrix<double> a1(n,1);
    Matrix<double> b1(n,1);
    Matrix<double> h(n,1);
    Matrix<double> xdraw(n,1);
    Matrix<int>    indx(n,1);
    Matrix<double> d_copy(n,n), dinv(n,n), tau(n,n), tau_(n,n), tau_copy(n,n), tauinv(n,n), c(n,n);

    d_copy = d;
    dinv = d_copy.inv();
    anu  = d*amu;

    tau_ = d*sigma;
    tau =  tau_*d.tran();

    tau_copy = tau;
    tauinv   = tau_copy.inv();

    a1 = a - anu;
    b1 = b - anu;
    for (int i=1; i<=n; i++)
    {
        aa = tauinv.get(i,i);
        h.set(i,1,1/std::sqrt(aa));
        for (int j=1; j<=n; j++)
            c.set(i,j, -tauinv.get(i,j)/aa);
    }

    for (int initer=0; initer<niter; initer++)
    {
        for (int k=1; k<=n; k++)
        {
            int i = kstep.get(k,1);
            double aa = 0;
            for (int j=1; j<=n; j++)
            {
                if (i != j)
                    aa = aa + c.get(i,j) * z.get(j,1);
            }

            if (la.get(i,1)==1)
                t1=normrt_rnd(0,1,(b1.get(i,1)-aa)/h.get(i,1));
            else if (lb.get(i,1)==1)
                t1=normlt_rnd(0,1,(a1.get(i,1)-aa)/h.get(i,1));
            else
                t1=normt_rnd(0,1,(a1.get(i,1)-aa)/h.get(i,1),(b1.get(i,1)-aa)/h.get(i,1));
            z.set(i,1,aa + h.get(i,1)*t1);
        }
    }
    //Transform back to x
    xdraw = dinv*z;
    for (int i=1; i <=n; i++)
    {
        xdraw.set(i,1,(xdraw.get(i,1) + amu.get(i,1)));
    }
    return xdraw;
}


/**** The mulit-variable Truncated Gaussian Distribution ****/
Matrix<double> Random::MvNormal(Matrix<double>& mean, Matrix<double>& covariance, Matrix<double>& MinValue, Matrix<double> & MaxValue)
{
    Matrix<double> vector(mean.row(), 1);
    Matrix<double> dec = covariance.chold();
    double rn_value;
    double left_value,right_value;

    int i;
    for (i=1; i<=(int)mean.row(); i++)
    {
        left_value  =  MinValue.get(i,1);
        right_value =  MaxValue.get(i,1);

        rn_value = normt_rnd(0,1,left_value,right_value);
        vector.set(i, 1, rn_value);
    }
    vector = dec*vector + mean;

    return vector;
}

Matrix<double> Random::MvNormal(Matrix<double>& mean, Matrix<double>& covariance, int SampleNum)
{
    Matrix<double> vector(mean.row(), SampleNum);

    for(int i=1; i<=SampleNum; i++)
    {
        vector.setColV(i,MvNormal(mean,covariance));
    }

    return vector;
}

Matrix<double> Random::MvNormal(Matrix<double>& mean, Matrix<double>& covariance, Matrix<double>& MinValue, Matrix<double> & MaxValue, int SampleNum)
{
    Matrix<double> vector(mean.row(), SampleNum);

    for(int i=1; i<=SampleNum; i++)
    {
        vector.setColV(i,MvNormal(mean,covariance,MinValue,MaxValue));
    }

    return vector;
}
double Random::randn()
{
    return gaussian(0.0, 1.0);
}

/**** The univariate Skew Gaussian Distribution ****/
double Random::rsn(double location, double scale, double shape)
{
    double u1 = randn();
    double u2 = randn();
    if (u2 > shape * u1)
    {
        u1 = -u1;
    }
    double r = location + scale * u1;
    return r;
}
double Random::dsn(double x, double location, double scale, double shape)
{
    double d = 2 * gaussian_pdf((x - location) / scale, 0.0, 1.0)
               * normal_01_cdf((shape * (x - location) / scale)) / scale;
    return d;
}
double Random::randu()
{
    return genrand_real3();
}
/**** The Gaussian Distribution PDF****/
double Random::gaussian_pdf(const double x, const double mu, const double sigma)
{
    double u;
    double p;

    u = (x - mu) / abs(sigma);
    p = (1 / (std::sqrt(2 * PI) * abs(sigma))) * std::exp(-u * u / 2);

    return p;
}
}
