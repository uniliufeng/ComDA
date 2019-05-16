#ifndef __LDAS_UnscentedKalmanFilter_H
#define __LDAS_UnscentedKalmanFilter_H

#include "filter.h"

using namespace arma;
namespace ldas
{
class UnscentedKalmanFilter : public Filter
{
public:
    /** Default constructor */
    UnscentedKalmanFilter();
    /** Default destructor */
    virtual ~UnscentedKalmanFilter();
    /** Copy constructor
     *  \param other Object to copy from
     */
    UnscentedKalmanFilter(const UnscentedKalmanFilter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    UnscentedKalmanFilter& operator=(const UnscentedKalmanFilter& other);

    /** constructor with parameters
    * \param state, unsigned int, the state number
    * \param obs, unsigned int, the observation number
    * \param a, double, default is 1.0, alpha parameter
    * \param b, double, default is 2.0, beta parameter
    * \param k, double, default is 0, kappa parameter
    */
    UnscentedKalmanFilter(const unsigned int state, const unsigned int obs, const double a=1.0, const double b=2.0, const double k=0);

    ///get state estimation
    mat getXa() const;
    /** Sigma points around reference point
    * \param x, mat(state_num,1), the reference point
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point
    * \param X, mat(state_num,state_num*2+1), output matrix with sigma conversion
    */
    void sigmas(const mat& x,const mat& P,  mat& X);
    /** Unscented Transformation
    * \param Y, mat, sigma points
    * \param Wm, mat, weights for mean
    * \param Wc, mat, weights for covraiance
    * \param y, mat, output matrix, transformed mean
    * \param Y1, mat, output matrix, transformed deviations
    */
    void ut(const mat& Y,const mat& Wm,const mat& Wc, mat& y, mat& Y1);
    /** update with observation
    * \param Xf, mat(state_num,sigma_num), the forecast state matrix
    * \param Yf, mat(observe_num,sigma_num), the forecast observation matrix
    * \param Yo, mat(observe_num,1), the observation matrix
    * \param Q, mat(state_num,state_num), the model error additive coveriance matrix
    * \param R, mat(observe_num,observer_num), the error additive coveriance matrix for observation
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point, output matrix
    */
    void update(const mat& Xf,const mat& Yf, const mat& Yo, const mat& Q, const double R, mat& P);
    /** update without observation
    * \param Xf, mat(state_num,sigma_num), the forecast state matrix
    * \param Q, mat(state_num,state_num), the model error additive coveriance matrix
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point
    */
    void update(const mat& Xf, const mat& Q, mat& P);
    /// set alpha parameter
    double alpha() const
    {
        return m_alpha;
    }
    /// get alpha parameter
    void alpha(const double a)
    {
        m_alpha=a;
    }
    /// set beta parameter
    double beta() const
    {
        return m_beta;
    }
    /// get beta parameter
    void beta(const double b)
    {
        m_beta=b;
    }
    /// set kappa parameter
    double kappa() const
    {
        return m_kappa;
    }
    /// get kappa parameter
    void kappa(const double k)
    {
        m_kappa=k;
    }
    /** compute RMSE
    * \param xtrue, mat(state_num,1), true value of the state matrix
    * \return double, the RMSE matrix for each state
    */
    double rmse(const mat& xtrue) const;
protected:
    double m_alpha,m_beta,m_kappa;
    mat Xa;
private:
};
typedef UnscentedKalmanFilter UKF;
}
#endif // __LDAS_UnscentedKalmanFilter_H
