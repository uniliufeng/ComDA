#ifndef __LDAS_EKF_H
#define __LDAS_EKF_H

#include "filter_itpp.h"

namespace ldas
{
class EKF : public Filter_itpp
{
public:
    /** Default constructor */
    EKF();
    /** Default destructor */
    virtual ~EKF();
    /** Copy constructor
     *  \param other Object to copy from
     */
    EKF(const EKF& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    EKF& operator=(const EKF& other);

	/** \brief 构造函数
	 * \param statedim 状态数
	 * \param obsdim 观测数
	 * \param EnNum 状态变量
	 * \param m_M 模型算子
	 * \param m_H 观测算子
	 */
	EKF(const itpp::vec xf_init, const itpp::vec obs, const int xdim, const int odim, itpp::mat C_, itpp::mat Sstate_Cov, itpp::mat Model_Cov, itpp::mat Obs_Cov);

	///更新
	void update();
	itpp::vec GetXa() const;
protected:
private:
    int Xdim; // extract state dimension
    int Odim; // extract observation dimension
    int U2dim; // extract exogenous input 2 dimension
    int Vdim; // extract process noise dimension
    int Ndim; // extract observation noise dimension
    itpp::vec xf; // copy particle buffer
    itpp::vec Obs;
    itpp::vec Hx;
    itpp::mat Px_; // model state forecast (Xdim*1)
    itpp::mat C;	//C = dh/dx
    itpp::vec Xa; // model state analysis (Xdim*1)
    itpp::vec Model_Noise_Mu;
    itpp::vec Obs_Noise_Mu;
    itpp::mat Sstate;
    itpp::mat Model_Covariance;
    itpp::mat Obs_Covariance;
};
}
#endif // __LDAS_EKF_H
