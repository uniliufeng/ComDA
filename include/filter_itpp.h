#ifndef __LDAS_FILTER_ITPP_H
#define __LDAS_FILTER_ITPP_H

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/itsignal.h>
#include "constant.h"
namespace ldas
{
class Filter_itpp
{
public:
    /** Default constructor */
    Filter_itpp();
    /** Default destructor */
    virtual ~Filter_itpp();
    /**
    GAUSEVAL  Calculates the likelihood of a dataset X under a Gaussian probability
              density described by the Gaussian data structure 'gausDS',
              i.e. P(X|gausDS). The column vectors of X are treated as IID samples to be
              evaluated. The function return a likelihood row-vector that contains one
              likelihood value for each IID sample in X.

      likelihood = gauseval(gausDS, X)

      INPUT
             gausDS       Gaussian data structure with the following fields
                .cov_type (string) covariance matrix type 'full' , 'diag' , 'sqrt' , 'sqrt-diag'
                .dim      (scalar)   dimension
                .mu       (c-vector) mean vector  (dim-by-1)
                .cov      (matrix)   covariance matrix of type cov_type  (dim-by-dim)
             X            (dim-by-M) buffer of M dim-by-1 data set vectors to be evaluated
             logflag      <optional> if 'logflag==1' then the probability is calculated in the log
                                     domain, i.e. likelihood=log(likelihood) (default : logflag=0)
      OUTPUT
             likelihood   (1-by-M) buffer of likelihood values
    */
    itpp::vec gauseval(int dim, itpp::vec mu, itpp::mat &cov, itpp::mat &X);
    ///element-by-element powers for vector
    itpp::vec power(const itpp::vec &v1,const itpp::vec &v2);
    itpp::vec likelihood(itpp::mat &obs, itpp::mat &state);
    itpp::ivec residualResample(const itpp::ivec &inIndex,const itpp::vec &weights);
    itpp::ivec fix(const itpp::vec &v);
    itpp::vec cumprod(const itpp::vec &v);
    void qr_economy(const itpp::mat &A, itpp::mat &Q, itpp::mat &R);
    itpp::mat cholupdate(const itpp::mat &R, const itpp::vec &X_, const char &ch);
protected:
private:
	void rotg(double& a,double& b,double& c,double& s);
	void dchud(int p, const itpp::vec &x, itpp::vec &c, itpp::vec &s, itpp::mat &r);
	double dot_r(const int N, const double *X, const int incX, const double *Y, const int incY);
	double dnrm2( const int N, const double *X, const int incX);
	int dchdd(int p, const itpp::vec &x, itpp::vec &c, itpp::vec &s, itpp::mat &r);

};
}
#endif // __LDAS_FILTER_ITPP_H
