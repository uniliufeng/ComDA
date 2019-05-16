#include "rng.h"

using namespace ldas;
Rng::Rng():Random()
{
    //ctor
}

Rng::~Rng()
{
    //dtor
}

itpp::vec Rng::randn(int size)
{
    itpp::vec temp(size);
    for (int i = 0; i < size; i++)
    {
        temp(i) = Random::randn();
    }
    return temp;
}
itpp::mat Rng::randn(int rows, int cols)
{
    itpp::mat temp(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            temp(i, j) = Random::randn();
        }
    }
    return temp;
}

/**** The Multivariate Skew Gaussian Distribution ****/
itpp::mat Rng::rmsn(int n, itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha)
{
    int k = Omega.cols();
    sn_struct Z = msn_quantities(xi, Omega, alpha);
    itpp::mat y = randn(n, k) * chol(Z.Psi);
    //each row of y is N_k(0,Psi)
    itpp::vec abs_y0_vec = abs(randn(n));
    itpp::mat abs_y0 = repeat(reshape(abs_y0_vec, n, 1), k);
    itpp::vec delta = Z.delta;
    itpp::mat z = elem_mult(repeat(reshape(delta, k, 1), n), transpose(abs_y0))
                  + elem_mult(repeat(reshape(sqrt(1 - pow(delta, 2)), k, 1), n),
                              transpose(y));
    y = transpose(repeat(reshape(xi, k, 1), n) + elem_mult(repeat(reshape(
                      (Z.omega_rebel), k, 1), n), z));
    return transpose(y);
}
itpp::vec Rng::dmsn(itpp::mat &x, itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha)
{
    itpp::vec scale = sqrt(diag(Omega));
    int k = x.rows();
    int n = x.cols();
    itpp::mat X = x - repeat(reshape(xi, length(xi), 1), n);
    itpp::mat z = elem_div(X, repeat(reshape(scale, length(scale), 1), n));
    itpp::vec Q = diag(transpose(X) * inv(Omega) * X); //diagonal of (x Omega^(-1) x^T)
    double Det = det(Omega);
    itpp::vec temp = transpose(z) * alpha;
    itpp::vec pdf = elem_mult(2 * exp(-0.5 * Q), normcdf(temp)) / std::sqrt(pow((2
                    * 3.1415926535897), k) * Det);
    return pdf;
}

itpp::vec Rng::uniform_vec(const double a, const double b, int size)
{
    itpp::vec temp(size);

    for (int i = 0; i < size; i++)
    {
        temp(i) = uniform(a, b);
    }
    return temp;
}

itpp::vec Rng::randu(int size)
{
    itpp::vec temp(size);
    for (int i = 0; i < size; i++)
    {
        temp(i) = Random::randu();
    }
    return temp;
}
itpp::mat Rng::randu(int rows, int cols)
{
    itpp::mat temp(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            temp(i, j) = Random::randu();
        }
    }
    return temp;
}

itpp::mat Rng::laplace_mat(itpp::vec &mu, itpp::vec &b, int rows, int cols)
{
    itpp::mat temp(rows, cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            temp(i, j) = laplace(b(i)) + mu(i);
        }
    }
    return temp;
}

/**** The Multivariate Gaussian Distribution ****/
/*Random vectors from a multivariate Normal distribution
 *CALL:  r = mnormrnd(n,mu,S)
 * r = matrix of random numbers from the multivariate normal
 *     distribution with mean mu and covariance matrix S.
 * n = number of sample vectors
 * method = 'svd'  Singular value decomp.  (stable, quite fast)
 * S must be a symmetric, semi-positive definite matrix with size equal
 * to the length of mu.
 * Example
 * mu = "0, 5"; S = "1 0.45; 0.45 0.25";
 * r = mnormrnd(100,mu,S)
 */
itpp::mat Rng::mnormrnd(itpp::vec &mu, itpp::mat &sigma, int cases)
{
    int m, n;
    m = sigma.rows();
    n = sigma.cols();
    if (m != n)
        cerr << "sigma must be square" << endl;

    int rows = mu.size();
    if (m != rows)
        cerr << "The length of mu must equal the number of rows in sigma"
             << endl;

    itpp::vec S_vec;
    itpp::mat U, V, S;
    svd(sigma, U, S_vec, V);

    S = diag(S_vec, 0);
    int U_row, V_row;
    U_row = U.rows();
    V_row = V.rows();

    int L = n;

    itpp::mat T = U(0, U_row - 1, 0, L - 1) * sqrt(S(0, L - 1, 0, L - 1))
                  * transpose(V(0, V_row - 1, 0, L - 1)); //squareroot of S


    itpp::mat r(m, cases), mu_matrix(m, cases);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < cases; j++)
        {
            mu_matrix(i, j) = mu(i);
        }
    }
    //cout <<mu_matrix<<endl;
    //cout <<prod(T,randn)<<endl;
    r = T * randn(m, cases) + mu_matrix;
    //cout<<r<<endl;
    return r;
}

/**** The Truncated Multivariate Gaussian Distribution ****/
/*Procedure for drawing from truncated multivariate normal based on
 Geweke's code. i.e draws from
 xdraw is N(amu,sigma) subject to a < x < b
 where N(.,) in the n-variate Normal, a and b are nx1
 la and lb are nx1 vectors set to one if no upper/lower bounds
 C  la(N)        If .TRUE., no corresponding lower bound
 C  lb(N)        If .TRUE., no corresponding upper bound
 */
itpp::vec Rng::tnorm_rnd(itpp::vec &amu, itpp::mat &sigma, itpp::vec &a, itpp::vec &b, itpp::vec &la, itpp::vec &lb)
{
    int niter = 10;
    int n = amu.size();
    itpp::mat d = itpp::eye(n); //Matrix of linear combinations of X for constraints:subject to a < d*x < b

    itpp::ivec kstep(n);
    for (int i = 0; i < n; i++)
    {
        kstep(i) = i; //order of Gibbs within the constraint rows
    }
    //transform to work in terms of z=d*x
    itpp::vec z = itpp::zeros(n);
    itpp::mat dinv = inv(d);
    itpp::vec anu = d * amu;

    itpp::mat tau = d * sigma * transpose(d);
    itpp::mat tauinv = inv(tau);
    itpp::vec a1 = a - anu;
    itpp::vec b1 = b - anu;
    itpp::mat c = itpp::zeros(n, n);
    itpp::vec h = itpp::zeros(n);

    double aa, t1;

    for (int i = 0; i < n; i++)
    {
        aa = tauinv(i, i);
        h(i) = 1.0 / std::sqrt(aa);
        for (int j = 0; j < n; j++)
        {
            c(i, j) = -tauinv(i, j) / aa;
        }
    }

    for (int initer = 0; initer < niter; initer++)
    {
        for (int i1 = 0; i1 < n; i1++)
        {
            int i = kstep(i1);
            aa = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    ;
                {
                    aa = aa + c(i, j) * z(j);
                }
            }

            if (la(i) == 1)
            {
                t1 = normrt_rnd(0, 1, (b1(i) - aa) / h(i));
            }
            else if (lb(i) == 1)
            {
                t1 = normlt_rnd(0, 1, (a1(i) - aa) / h(i));
            }
            else
            {
                t1 = normt_rnd(0, 1, (a1(i) - aa) / h(i), (b1(i) - aa) / h(i));
            }

            z(i) = aa + h(i) * t1;
        }
    }

    //Transform back to x
    itpp::vec xdraw = dinv * z;
    for (int i = 0; i < n; i++)
    {
        xdraw(i) = xdraw(i) + amu(i);
    }

    return xdraw;
}

/**** The latin hypercube sample of the Multivariate Gaussian Distribution ****/
/*%LHSNORM Generate a latin hypercube sample with a normal distribution
 %   X=LHSNORM(MU,SIGMA,N) generates a latin hypercube sample X of size
 %   N from the multivariate normal distribution with mean vector M
 %   and covariance matrix SIGMA.  X is similar to a random sample from
 %   the multivariate normal distribution, but the marginal distribution
 %   of each column is adjusted so that its its sample marginal
 %   distribution is close to its theoretical normal distribution.
 * */
itpp::mat Rng::lhsnorm(itpp::vec &mu, itpp::mat &sigma, int n)
{
    // Generate a random sample with a specified distribution and
    // correlation structure -- in this case multivariate normal
    itpp::mat z = transpose(mnormrnd(mu, sigma, n));

    // Find the ranks of each column
    int p = length(mu);
    itpp::mat x = itpp::zeros(z.rows(), z.cols());
    itpp::vec x_col, z_col;
    for (int i = 0; i < p; i++)
    {
        z_col = z.get_col(i);
        x_col = rank(z_col);
        x.set_col(i, x_col);
    }

    // Get gridded or smoothed-out values on the unit interval
    x = x - randu(x.rows(), x.cols());

    x = x / n;

    // Transform each column back to the desired marginal distribution,
    // maintaining the ranks (and therefore rank correlations) from the
    // original random sample
    for (int i = 0; i < p; i++)
    {
        x_col = x.get_col(i);
        x_col = norminv(x_col, mu(i), std::sqrt(sigma(i, i)));
        x.set_col(i, x_col);
    }

    return transpose(x);
}
itpp::mat Rng::lhsnorm_trunctated(itpp::vec &mu, itpp::mat &sigma, int n, itpp::vec &a, itpp::vec &b, itpp::vec &la,
                                  itpp::vec &lb)
{
    // Generate a random sample with a specified distribution and
    // correlation structure -- in this case multivariate normal
    itpp::vec z_row;
    itpp::mat z(n, length(mu));
    for (int i = 0; i < n; i++)
    {
        z_row = tnorm_rnd(mu, sigma, a, b, la, lb);
        z.set_row(i, z_row);
    }

    // Find the ranks of each column
    int p = length(mu);
    itpp::mat x = itpp::zeros(z.rows(), z.cols());
    itpp::vec x_col, z_col;
    for (int i = 0; i < p; i++)
    {
        z_col = z.get_col(i);
        x_col = rank(z_col);
        x.set_col(i, x_col);
    }

    // Get gridded or smoothed-out values on the unit interval
    x = x - randu(x.rows(), x.cols());

    x = x / n;

    // Transform each column back to the desired marginal distribution,
    // maintaining the ranks (and therefore rank correlations) from the
    // original random sample
    for (int i = 0; i < p; i++)
    {
        x_col = x.get_col(i);
        x_col = norminv(x_col, mu(i), std::sqrt(sigma(i, i)));
        x.set_col(i, x_col);
    }

    return transpose(x);
}

itpp::vec Rng::normcdf(itpp::vec &x)
{
    itpp::vec temp(length(x));
    for (int i = 0; i < x.size(); i++)
    {
        temp(i) = normal_01_cdf(x(i));
    }
    return temp;
}

itpp::vec Rng::rank(itpp::vec &x)
{
    // Similar to tiedrank, but no adjustment for ties here
    itpp::ivec rowidx = sort_index(x);
    itpp::vec r(length(x));
    for (int i = 1; i <= length(x); i++)
    {
        r(rowidx(i - 1)) = i;
    }

    return r;
}

itpp::vec Rng::norminv(itpp::vec &p, double mu, double sigma)
{
    itpp::vec p_ = 2 * p;

    itpp::vec x0 = -std::sqrt(2) * erfcinv(p_);

    return sigma * x0 + mu;
}

itpp::vec Rng::erfcinv(itpp::vec &y)
{
    return erfinv(1 - y);
}

sn_struct Rng::msn_quantities(itpp::vec &xi, itpp::mat &Omega, itpp::vec &alpha)
{
    sn_struct output;

    int k = length(alpha);

    itpp::vec omega = sqrt(diag(Omega));
    itpp::mat O_cor = diag(1.0 / omega) * Omega * diag(1.0 / omega);
    double tmp = std::sqrt(1 + alpha * (O_cor * alpha));
    itpp::vec delta = O_cor * alpha / tmp;
    itpp::vec lambda = elem_div(delta, sqrt(1 - pow(delta, 2)));
    itpp::mat D = diag(sqrt(1 + pow(lambda, 2)));
    itpp::mat Psi = D * (O_cor - (reshape(delta, k, 1) * reshape(delta, 1, k))) * D;
    Psi = (Psi + transpose(Psi)) / 2.0;
    itpp::mat O_inv = inv(Omega);
    itpp::vec oi = sqrt(diag(O_inv));
    itpp::mat O_conc = diag(1.0 / oi) * (-O_inv) * diag(1.0 / oi);
    O_conc = O_conc - diag(diag(O_conc)) + itpp::eye(k);//1 on the dagonal of O_conc
    itpp::vec muZ = delta * std::sqrt(2 / 3.1415926535897);
    itpp::vec muY = xi + elem_mult(omega, muZ);
    itpp::mat Sigma = diag(omega) * (O_cor - reshape(muZ, k, 1) * reshape(muZ, 1, k))
                      * diag(omega);
    Sigma = (Sigma + transpose(Sigma)) / 2.0;
    itpp::vec cv = elem_div(muZ, sqrt(1 - pow(muZ, 2)));
    itpp::vec gamma1 = 0.5 * (4 - 3.1415926535897) * pow(cv, 3);

    output.xi = xi;
    output.Omega = Omega;
    output.alpha = alpha;
    output.omega_rebel = omega;
    output.mean = muY;
    output.variance = Sigma;
    output.Omega_cor = O_cor;
    output.Omega_conc = O_conc;
    output.Psi = Psi;
    output.lambda = lambda;
    output.delta = delta;
    output.skewness = gamma1;

    return output;
}

