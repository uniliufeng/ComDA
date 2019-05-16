#ifndef __LDAS_SPPF_H
#define __LDAS_SPPF_H

#include "filter_itpp.h"
#include "constant.h"

namespace ldas
{
	///SPPF  Sigma-Point Particle Filter.
class SPPF : public Filter_itpp
{
public:
    /** Default constructor */
    SPPF();
    /** Default destructor */
    virtual ~SPPF();
    /** Copy constructor
     *  \param other Object to copy from
     */
    SPPF(const SPPF& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    SPPF& operator=(const SPPF& other);


	/** \brief 构造函数
	 * \param statedim 状态数
	 * \param obsdim 观测数
	 * \param EnNum 状态变量
	 * \param m_M 模型算子
	 * \param m_H 观测算子
	 */
	SPPF(const itpp::mat particles_init, itpp::mat particlesNew_init, itpp::mat particlesPred_init, itpp::vec prior_init, itpp::Array < itpp::mat > SxPred_init, const itpp::vec obs, const float resampleThreshold);

	///更新
	void update();
	itpp::vec GetXa() const;
	itpp::mat GetXaEn() const;

protected:
private:
    int Xdim; // extract state dimension
    int Odim; // extract observation dimension
    int U2dim; // extract exogenous input 2 dimension
    int Vdim; // extract process noise dimension
    int Ndim; // extract observation noise dimension
    int N; // number of particles
    itpp::mat particles; // copy particle buffer
    itpp::mat particlesNew;
    itpp::mat particlesPred;
    itpp::Array < itpp::mat > SxPred;
    itpp::vec Obs;
    itpp::vec weights; // particle weights
    itpp::vec prior; // calculate transition prior for each particle (in log domain)
    int St; // resample threshold
    itpp::vec Xa; // mean of model state analysis (Xdim*1)
};
}
#endif // __LDAS_SPPF_H
