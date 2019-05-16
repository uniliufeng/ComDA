#ifndef __LDAS_UNSCENTEDPARTICLEFILTER_H
#define __LDAS_UNSCENTEDPARTICLEFILTER_H

#include "particlefilter.h"

namespace ldas
{
class UnscentedParticleFilter : public ParticleFilter
{
public:
    /** Default constructor */
    UnscentedParticleFilter();
    /** Default destructor */
    virtual ~UnscentedParticleFilter();
    /** Copy constructor
     *  \param other Object to copy from
     */
    UnscentedParticleFilter(const UnscentedParticleFilter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    UnscentedParticleFilter& operator=(const UnscentedParticleFilter& other);
    /** constructor with parameters
    * \param StateNum, unsigned int, the state number
    * \param ObsNum, unsigned int, the observation number
    * \param EnNum, unsigned int, the particles number
    * \param a, double, default is 1.0, alpha parameter
    * \param b, double, default is 2.0, beta parameter
    * \param k, double, default is 0, kappa parameter
    */
    UnscentedParticleFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum, const double a=1.0, const double b=2.0, const double k=0);
    /** update with observation
    * \param state_en_sigmas, mat(state_num,(state_num*2+1)*particle_num), the forecast state matrix, M(X)
    * \param obs_forecast_en_sigmas, mat(observe_num,(state_num*2+1)*particle_num), the forecast observation matrix, H(X)
    * \param obs, mat(observe_num,1), the observation matrix
    * \param R, mat(observe_num,observe_num), the error additive coveriance matrix for observation
    * \param P[(state_num*2+1)*particle_num], the coveriance matrix for each sigma point.
    */
    void update(const mat& state_en_sigmas, const mat& obs_forecast_en_sigmas, const mat& obs, const mat& R,mat* P);
    /** update without observation
    * \param state, mat(state_num,particle_num), the forecast state matrix
    * \return void
    */
    void update(const mat& state/*, mat* */);
    /** unscented kalman filter with observation
    * \param Xf, mat(state_num,sigma_num), the forecast state matrix
    * \param Yf, mat(observe_num,sigma_num), the forecast observation matrix
    * \param Yo, mat(observe_num,1), the observation matrix
    * \param R, mat(observe_num,observer_num), the error additive coveriance matrix for observation
    * \param Xa, mat(state_num,1), output state matrix
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point
    */
    void ukf(const mat& Xf, const mat& Yf, const mat& Yo, const mat& R,mat& Xa, mat& P);
    /** unscented kalman filter without observation
    * \param Xf, mat(state_num,sigma_num), the forecast state matrix
    * \param Xa, mat(state_num,1), output state matrix
    * \param P, mat(state_num,state_num), the coveriance matrix for each sigma point
    */
    void ukf(const mat& Xf, mat& Xa, mat& P);
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
protected:
    double m_alpha,m_beta,m_kappa;
private:
};
typedef UnscentedParticleFilter UPF;
}
#endif // __LDAS_UNSCENTEDPARTICLEFILTER_H
