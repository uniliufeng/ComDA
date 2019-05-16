#ifndef __LDAS_PARTICLEFILTER_H
#define __LDAS_PARTICLEFILTER_H

#include "constant.h"
#include "filter.h"

namespace ldas
{
///General Particle Filter, including the SIS and SIR
class ParticleFilter: public Filter
{
public:
    /** Default constructor */
    ParticleFilter();

    /** constructor with state, observation, particles number
     * \param StateNum, unsigned int, the state variable number
     * \param ObserveNum, unsigned int, the observe variable number
     * \param EnNum, unsigned int, the particles number
     */
    ParticleFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum);

    /** Default destructor */
    virtual ~ParticleFilter();

    /** update with observation
    * \param state, mat(state_num,particle_num), the forecast state particles
    * \param Yf, mat(observe_num,particle_num), the forecast observation matrix
    * \param obs, mat(observe_num,1), the observe matrix
    */
    void update(const mat& state, const mat& Yf,const mat& obs);
    /** update without observation
    * \param state, mat(state_num,particle_num), the forecast state particles
    */
    void update(const mat& state);
    ///get state estimation
    mat getXa() const;
    ///get state particles
    mat getXaEn() const;

    /** Observation likelihood function
    calculates p(y|x) for a given realization of the state variable 'state' and a particular observation instance 'obs'.
    i.e. Calculates the value of P(OBS|STATE) = P(y|x)
    * \param	obs, mat(observe_num*particle_num), observation at time k
    * \param	obs_forecast, mat(state_num*particle_num), the observe forecast from H function at time k
    * \return mat(1*particle_num), p(y(k)|x(k))
    */
    mat likelihood(const mat &obs, const mat &obs_forecast);

    /** set particles number
    * \param n, unsigned int, particles number
    */
    void particleNum(const unsigned int n);
    /// get particles number
    unsigned int particleNum() const;
    /// set resample threshold
    void threshold(const double th);
    /// get resample threshold
    double threshold() const;

    /** resample the states particles with weights
    * \param weights, mat(1*particle_num)
    * \param particles, old particles, mat(state_num*particle_num)
    * \return resampled particles, mat(state_num*particle_num)
    */
    mat residual_resample(const mat& weights,const mat& particles);
    /** compute gaussian likelihood pdf
    * \param dim, int, gaussian noise dimension or observe_num?
    * \param mu, mat, mean(expected? would be zero default)
    * \param cov, mat, coveriance
    * \param X, mat(state_num,particle_num)
    * \return mat(1,particle_num)
    */
    mat gauseval(const int dim, const mat& mu, const mat& cov, const mat& X);

protected:
    unsigned int particle_num;
    /// resample threshold
    double m_threshold;
    mat weights;
    mat particles; // copy particle buffer
    mat Xa; // mean of model state analysis (Xdim*1)

    ///initialize the particle filter
    void init();
private:
};
typedef ParticleFilter PF;
}
#endif
