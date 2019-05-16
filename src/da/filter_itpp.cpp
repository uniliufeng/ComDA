#include "filter_itpp.h"
#include "rng.h"

using namespace ldas;
Filter_itpp::Filter_itpp()
{
    //ctor
}

Filter_itpp::~Filter_itpp()
{
    //dtor
}

itpp::vec Filter_itpp::gauseval(int dim, itpp::vec mu, itpp::mat &cov, itpp::mat &X)
{
    double normfact;
    int row = X.rows();
    int col = X.cols(); // number of input vectors
    itpp::vec likelihood = itpp::zeros(col); // preallocate likelihood matrix
    itpp::mat XX, S, foo;
    // calculations depend on covariance type
    normfact = std::pow(2 *PI_NUMBER, dim / 2.0);
    XX = X - repeat(reshape(mu, row, 1), col);
    S = transpose(chol(cov));
    foo = inv(S) * XX;
    likelihood = exp(-0.5 * sum(elem_mult(foo, foo), 1)) / (normfact * abs(
                     prod(diag(S))));
    return likelihood;
}

itpp::vec Filter_itpp::power(const itpp::vec &v1,const itpp::vec &v2)
{
    itpp::vec out(v1.size());
    for(int i=0; i<v1.size(); i++)
    {
        out(i) = std::pow(v1(i),v2(i));
    }
    return out;
}

itpp::vec Filter_itpp::likelihood(itpp::mat &obs, itpp::mat &state)
/****************************************************************************************************
 % LIKELIHOOD  Observation likelihood function
 %
 % Function-handle to the observation likelihood function that calculates p(y|x) for a
 % given realization of the state variable 'state' and a particular observation instance 'obs'.
 %
 %   i.e. Calculates the value of P(OBS|STATE) = P(y|x)
 %
 %   INPUT
 %         model          GSSM data structure
 %         obs            observation at time k
 %         state          state at time k
 %         U2             exogeneous input to HFUN at time k
 %         oNoiseDS       (optional) measurement noise NoiseDS data structure to use for evaluation of
 %                        transition prior. If this is ommitted, model.oNoise, is used.
 %   OUTPUT
 %         llh            p(y(k)|x(k))
 %
 %-- This function must be defined by the user!
 ******************************************************************************************************/
{
    int dim = state.rows();
    int obsdim = obs.rows();
    int nov = state.cols();
    itpp::vec llh(nov);
    itpp::vec zero = itpp::zeros(obsdim);
    itpp::vec X_col(obsdim), state_col(dim), obs_col(obsdim);
    itpp::mat X(obsdim, nov);
    for (int i = 0; i < nov; i++)
    {
        obs_col = obs.get_col(i);
        state_col = state.get_col(i);
        itpp::vec obssim = state_col;
        if (obsdim > 1)
        {
            for (int j = 0; j < obsdim; j++)
            {
                X_col(j) = obs_col(j) - obssim(0);
            }
        }
        else
        {
            X_col = obs_col - obssim;
        }

        X.set_col(i, X_col);
    }
    itpp::vec mu = itpp::zeros(dim);
    mu = 0.0;
    itpp::mat covariance = itpp::zeros(dim,dim);
    for(int i=0; i<dim; i++)
    {
        covariance(i,i) = variance(state.get_row(i));
    }
    llh = gauseval(dim, mu, covariance, X);
    return llh;
}
itpp::ivec Filter_itpp::residualResample(const itpp::ivec &inIndex,const itpp::vec &weights)
{
    int S = length(weights);		// S = Number of particles

    itpp::ivec outIndex(S);	// setup output index buffer

    //=== RESIDUAL RESAMPLING  ==========================================================
    itpp::ivec N_kind(S);

    // first integer part
    itpp::vec weights_res = S * weights;
    N_kind = fix(weights_res);

    // residual number of particles to sample
    int N_res = S - sum(N_kind);
    Rng rng;

    if(N_res)
    {
//		cout<<N_res<<endl;
        weights_res = (weights_res - N_kind) / N_res;
        itpp::vec cumDist = cumsum(weights_res);
        // generate N_res ordered random variables uniformly distributed in [0,1]

        itpp::vec temp1(N_res);
        for(int i=0; i<N_res; i++)
        {
            temp1(i) = N_res - i;
        }

        itpp::vec randvec,powervec,divec,u;
        randvec = rng.randu(N_res);
        divec = 1.0/temp1;
        powervec = power(randvec,divec);
        u = reverse(cumprod(powervec));

        for(int i=0; i<N_res; i++)
        {
            int j = 0;
            while(u(i)>cumDist(j))
            {
                j = j+1;
            }
            N_kind(j) = N_kind(j) + 1;
        }
    }
//	cout<<N_kind<<endl;
    //=== COPY RESAMPLED TRAJECTORIES =====================================================
    int index = 0;
    for(int i=0; i<S; i++)
    {
        if(N_kind(i)>0)
        {
            int j;
            j = index;
            do
            {
                outIndex(j) = inIndex(i);
                j++;
            }
            while(j<=index+N_kind(i)-1);
        }
        index = index + N_kind(i);
    }
    return outIndex;
}
//rounds the elements of X to the nearest integers towards zero
itpp::ivec Filter_itpp::fix(const itpp::vec &v)
{
    itpp::ivec out(v.size());
    for(int i=0; i<v.size(); i++)
    {
        if(v(i)>0)
        {
            out(i) = int(floor(v(i)));
        }
        else
        {
            out(i) = int(ceil(v(i)));
        }
    }

    return out;
}

//Cumulative product of elements
itpp::vec Filter_itpp::cumprod(const itpp::vec &v)
{
    itpp::vec out(v.size());
    out(0)=v(0);
    for (int i=1; i<v.size(); i++)
    {
        out(i) = out(i-1) * v(i);
    }

    return out;
}

//where A is m-by-n,produces the "economy size" decomposition.
//If m>n, only the first n columns of Q and the first n rows of R are
//computed. If m<=n, this is the same as [Q,R] = QR(A).
void Filter_itpp::qr_economy(const itpp::mat &A, itpp::mat &Q, itpp::mat &R)
{
 /*   int m, n, k, info, lwork, i, j;

    m = A.rows();
    n = A.cols();
    lwork = 1000 * n;
    k = std::min(m, n);
    itpp::vec tau(k);
    itpp::vec work(lwork);

    if (m > n)
    {
        R = A;
        dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);
        Q = R;

        // construct R
        for (i = 0; i < m; i++)
            for (j = 0; j < std::min(i, n); j++)
                R(i, j) = 0;

        Q.set_size(m, m, true);
        dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
                &info);

        Q = Q.get_cols(0, n - 1);
        R = R.get_rows(0, n - 1);
    }
    else
    {
        R = A;
        dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);
        Q = R;

        // construct R
        for (i = 0; i < m; i++)
            for (j = 0; j < std::min(i, n); j++)
                R(i, j) = 0;

        Q.set_size(m, m, true);
        dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
                &info);
    }*/
}

void Filter_itpp::rotg(double& a,double& b,double& c,double& s)
{
    const double roe = (abs(a) > abs(b) ? a : b);
    const double scale = abs(a) + abs(b);
    double r, z;

    if (scale != 0.0)
    {
        const double aos = a / scale;
        const double bos = b / scale;
        r = scale * std::sqrt(aos * aos + bos * bos);
        r = (roe>=0?1:-1) * r;
        c = a / r;
        s = b / r;
        z = 1.0;
        if (abs(a) > abs(b))
            z = s;
        if (abs(b) >= abs(a) && c != 0.0)
            z = 1.0 / c;
    }
    else
    {
        c = 1.0;
        s = 0.0;
        r = 0.0;
        z = 0.0;
    }

    a = r;
    b = z;
}

//Perform a cholesky factor update
//Ported from the LINPACK FORTRAN function DCHUD
void Filter_itpp::dchud(int p, const itpp::vec &x, itpp::vec &c, itpp::vec &s, itpp::mat &r)
{
    int j, i, jm1;
    double t, xj, rtmp, ctmp, stmp;
    for (j = 1; j <= p; j++)
    {
        xj = x(j - 1);
        jm1 = j - 1;
        if (jm1 >= 1)
        {
            for (i = 1; i <= jm1; i++)
            {
                t = c(i - 1) * r(i - 1, j - 1) + s(i - 1) * xj;
                xj = c(i - 1) * xj - s(i - 1) * r(i - 1, j - 1);
                r(i - 1, j - 1) = t;
            }
        }
        rtmp = r(j - 1, j - 1);
        ctmp = c(j - 1);
        stmp = s(j - 1);
        rotg(rtmp, xj, ctmp, stmp);
        r(j - 1, j - 1) = rtmp;
        c(j - 1) = ctmp;
        s(j - 1) = stmp;
    }
}

double Filter_itpp::dot_r(const int N, const double *X, const int incX, const double *Y, const int incY)
{
    double r = 0.0;
    int i;
    int ix = ((incX) > 0 ? 0 : ((N) - 1) * (-(incX)));//OFFSET(N, incX);
    int iy = ((incY) > 0 ? 0 : ((N) - 1) * (-(incY)));//OFFSET(N, incY);

    for (i = 0; i < N; i++)
    {
        r += X[ix] * Y[iy];
        ix += incX;
        iy += incY;
    }

    return r;
}

double Filter_itpp::dnrm2( const int N, const double *X, const int incX)
{
    double  scale = 0.0, ssq   = 1.0;
    int   i,   ix    = 0;

    if (N <= 0 || incX <= 0)
        return 0;
    else if (N == 1)
        return abs(X[0]);

    for (i = 0; i < N; i++)
    {
        const double x = X[ix];

        if (x != 0.0)
        {
            const double ax = fabs(x);

            if (scale < ax)
            {
                ssq   = 1.0 + ssq * (scale / ax) * (scale / ax);
                scale = ax;
            }
            else
            {
                ssq += (ax / scale) * (ax / scale);
            }
        }

        ix += incX;
    }

    return scale * std::sqrt(ssq);
}

int Filter_itpp::dchdd(int p, const itpp::vec &x, itpp::vec &c, itpp::vec &s, itpp::mat &r)
{
    int info; //ldr,ldz,nz;
    int i, ii, j, k;
    double alpha, norm, a; //azeta,dnrm2;
    double t, xx, scale, b; //ddot,zeta;
    double tempVar;
    double rvectemp[20];
    double sVectemp[20];
    double cVectemp[20];
    info = 0;
    sVectemp[0] = x(0) / r(0, 0);
    if (p >= 2)
    {
        for (j = 2; j <= p; j++)
        {
            for (k = 0; k < p; k++)
            {
                rvectemp[k] = r(k, j - 1);
            }
            sVectemp[j - 1] = x(j - 1) - dot_r(j - 1, rvectemp, 1,
                                                    sVectemp, 1);
            sVectemp[j - 1] = sVectemp[j - 1] / r(j - 1, j - 1);
        }
    }
    for (k = 0; k < p; k++)
    {
    }
    norm = dnrm2(p, sVectemp, 1);
    if (norm < 1.0)
    {
        alpha = std::sqrt(1.0 - norm * norm);
        for (ii = 1; ii <= p; ii++)
        {
            i = p - ii + 1;
            scale = alpha + abs(sVectemp[i - 1]);
            a = alpha / scale;
            b = sVectemp[i - 1] / scale;
            norm = std::sqrt(a * a + b * b);
            cVectemp[i - 1] = a / norm;
            sVectemp[i - 1] = b / norm;
            alpha = scale * norm;
        }
        for (j = 1; j <= p; j++)
        {
            xx = 0;
            for (ii = 1; ii <= j; ii++)
            {
                i = j - ii + 1;
                t = cVectemp[i - 1] * xx + sVectemp[i - 1] * r(i - 1, j - 1);
                tempVar = cVectemp[i - 1] * r(i - 1, j - 1) - sVectemp[i - 1]
                          * xx;
                r(i - 1, j - 1) = tempVar;
                xx = t;
            }
        }

    }
    else
        info = -1;
    for (k = 0; k < p; k++)
    {
        s(k) = sVectemp[k];
        c(k) = cVectemp[k];
    }
    return info;
}

//Rank 1 update to Cholesky factorization.
itpp::mat Filter_itpp::cholupdate(const itpp::mat &R, const itpp::vec &X_, const char &ch)
{
    int LDR = R.rows();
    int P = LDR;
    int LDZ = 0;
    int NZ = 0;
    itpp::mat A(R);
    itpp::vec X(X_);
    itpp::vec Z(P * NZ), Y(NZ), RHO(NZ), C(P), S(P);

    switch (ch)
    {
    case '+':
        dchud(LDR, X, C, S, A);
        break;
    case '-':
        int err;
        err = dchdd(LDR, X, C, S, A);
        if (err == 0)
            //				{
            //					cout<<"the entire downdating was successful."<<endl;
            //				}
            //				else if(err==-1)
            //				{
            //					cerr<<"if R could not be downdated.  In this case, all quantities are left unaltered."<<endl;
            //				}
            //				else
            //				{
            //					cerr<<"if some RHO could not be downdated.  The offending RHO's are set to -1."<<endl;
            //				}
            break;
    default:
        std::cerr << "char is not a effective characater" << std::endl;
        break;
    }

    return A;
}
