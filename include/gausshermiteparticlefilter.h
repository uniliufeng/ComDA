#ifndef __LDAS_GAUSSHERIMTEPARTICLEFILTER_H
#define __LDAS_GAUSSHERIMTEPARTICLEFILTER_H

#include "particlefilter.h"
#include "gausshermitefilter.h"

namespace ldas
{
class GaussHerimteParticleFilter : public ParticleFilter,public GaussHermiteFilter
{
public:
    /** Default constructor */
    GaussHerimteParticleFilter();
    /** Default destructor */
    virtual ~GaussHerimteParticleFilter();
    /** Copy constructor
     *  \param other Object to copy from
     */
    GaussHerimteParticleFilter(const GaussHerimteParticleFilter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    GaussHerimteParticleFilter& operator=(const GaussHerimteParticleFilter& other);
    /** constructor with parameters
    * \param StateNum, unsigned int, the state number
    * \param ObsNum, unsigned int, the observation number
    * \param EnNum, unsigned int, the particles number
    */
    GaussHerimteParticleFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum);
    /** update with observation
    * \param state_en_gauss, mat(state_num,gauss_points*particle_num), the forecast state matrix via guass points transformation, M(X)
    * \param obs_forecast_en_gauss, mat(observe_num, gauss_points*particle_num), the forecast observation matrix, H(X)
    * \param obs, mat(observe_num,1), the observation matrix
    * \param R, mat(observe_num,observe_num), the error additive coveriance matrix for observation
    * \param P[particle_num], the coveriance matrix of the model for each particle
    */
    void update(const mat& state_en_gauss, const mat& obs_forecast_en_gauss, const mat& obs, const mat& R,mat* P);
    /** predict step
    * \param state_en_gauss, mat(state_num,gauss_points*particle_num), the forecast state matrix
    * \param Q, mat(state_num,state_num), model error coveriance matrix
    * \param P[particle_num], the coveriance matrix of the model for each particle
    */
    void update(const mat& state_en_gauss, const mat& Q, mat* P);
    /** particle filter update with observation
    * \param state, mat(state_num,particle_num), the forecast state particles
    * \param obs_en, mat(observe_num,particle_num), the forecast observation matrix
    * \param obs, mat(observe_num,1), the observe matrix
    */
    void update(const mat& state,const mat& obs_en,const mat& obs)
    {
    	ParticleFilter::update(state,obs_en,obs);
    }

    ///get state estimation
    mat getXa()
    {
    	return ParticleFilter::getXa();
    }
protected:
private:
};
typedef GaussHerimteParticleFilter GHPF;
}
#endif // __LDAS_GAUSSHERIMTEPARTICLEFILTER_H
