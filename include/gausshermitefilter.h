#ifndef __LDAS_GAUSSHERMITEFILTER_H
#define __LDAS_GAUSSHERMITEFILTER_H

#include "filter.h"

namespace ldas
{
class GaussHermiteFilter : public Filter
{
public:
    /** Default constructor */
    GaussHermiteFilter();
    /** Default destructor */
    virtual ~GaussHermiteFilter();
    /** Copy constructor
     *  \param other Object to copy from
     */
    GaussHermiteFilter(const GaussHermiteFilter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    GaussHermiteFilter& operator=(const GaussHermiteFilter& other);

    /** constructor with parameters
    * \param state, unsigned int, the state number
    * \param obs, unsigned int, the observation number
    * \param n, int, default is 3, gauss order
    */
    GaussHermiteFilter(const unsigned int state, const unsigned int obs);

    ///get state estimation
    mat getXa() const;
    /** Gauss points around reference point
    * \param x, mat(state_num,1), the reference point
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point
    * \param X, mat(state_num,3), output matrix with gauss conversion
    * \param flag, int, 1: use sqrt(diag(P)) way to compute P=trans(S)*S, else use cholesky way.
    */
    void gausspoint(const mat& x,const mat& P,  mat& X, int flag=1);

    /** update with observation
    * \param Xt, mat(state_num,gauss_points_num), Gauss points
    * \param Xf, mat(state_num,1), the forecast state matrix at current time(X_k|k-1)
    * \param Yt, mat(observe_num,gauss_points_num), the forecast observation matrix with gauss points, H(Xt)
    * \param Yo, mat(observe_num,1), the observation matrix
    * \param Q, mat(state_num,state_num), the model error additive coveriance matrix
    * \param R, mat(observe_num,observer_num), the error additive coveriance matrix for observation
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point, output matrix
    */
    void update(const mat& Xt,const mat& Xf,const mat& Yt,const mat& Yo,const mat& R,mat& P);
    /** update without observation
    * \param Xf, mat(state_num,sigma_num), the forecast state matrix
    * \param Q, mat(state_num,state_num), the model error additive coveriance matrix
    */
    void update(const mat& Xf, const mat& Q, mat& P);
    /** compute RMSE
    * \param xtrue, mat(state_num,1), true value of the state matrix
    * \return double, the RMSE matrix for each state
    */
    double rmse(const mat& xtrue) const;
protected:
	///state analysis result, dim: state_num*1
    mat Xa;
    ///model error variance, dim: state_num*state_num
    mat P;
    /// gauss order number, determine the result polynomial of degree up to 2*guass_order-1
    unsigned int gauss_order;
    /// gauss points number, = pow(gauss_order,state_num)
    unsigned int gauss_points_num;
    /// q parameter in the gauss points transformation, dim:gauss_order*1
    mat yi;
    /// q paramter for multi dimension state, dim:gauss_points_num*state_num
    mat qi;
    /// weights of the Gauss-Hermite quadrature, dim:gauss_order*1
    mat ai;
    /// weights for multi dimension state, equal to wi1*wi2*...*win, dim:gauss_points_num*state_num
    mat wi;
    /// wi production by row, pwi=prod(wi,1), dim: gauss_points_num*1
    mat pwi;
    /// initialize the ai,yi mat
    void init();
private:
};
typedef GaussHermiteFilter GHF;
}
#endif // __LDAS_GAUSSHERMITEFILTER_H
