#ifndef __LDAS_ENSEMBLEKALMANFILTER_H
#define __LDAS_ENSEMBLEKALMANFILTER_H

#include "filter.h"
#include "random.hpp"

using namespace arma;
namespace ldas
{
class EnsembleKalmanFilter : public Filter
{
public:
    /** Default constructor */
    EnsembleKalmanFilter();
    /** Default destructor */
    virtual ~EnsembleKalmanFilter();
    /** \brief 构造函数
    * \param StateNum 状态数
    * \param ObsNum 观测数
    * \param EnNum 集合数
    */
    EnsembleKalmanFilter(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum);

    /** Copy constructor
     *  \param other Object to copy from
     */
    EnsembleKalmanFilter(const EnsembleKalmanFilter& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    EnsembleKalmanFilter& operator=(const EnsembleKalmanFilter& other);
    /** update with observation.
    * \param XfEn n*N matrix, model forecast state ensemble
    * \param YEn m*N matrix, ensemble of estimated observation, with H(observe model)
    * \param Yo m*1 matrix, observation
    * \param R m*m, observe error coveriance
    * \param Q n*n, model error coveriance
    */
    void update(mat& XfEn, mat& YEn,mat& Yo,mat& R,mat& Q);
    void update(mat& XfEn);
    mat ensembleInc(const mat& ensemble);
    ///预报, generate XfEn, perturb XfEn.
    void predict(const mat& m_Q);
    void perturbYo(const mat& m_R);
    ///mean of model state analysis (n*1)
    mat GetXa() const;
    ///Matrix holding the Ensemble members forecast (n*N)
    mat GetXaEn() const;
    void is_perturb_model(const bool status);
    bool is_perturb_model() const;
    void is_perturb_observe(const bool status);
    bool is_perturb_observe() const;
protected:
    unsigned int ensemble_num;
private:
    /// mean of model state forecast (n*1)
    mat	Xf;
    /// mean of model state analysis (n*1)
    mat	Xa;
    /// Matrix holding the Ensemble members forecast (n*N)
    mat	XfEn;
    /// Matrix holding the Ensemble members analysis (n*N)
    mat  XaEn;
    /// mean of estimated observation (m*1)
    //mat	Y;
    /// ensemble of estimated observation (m*N)
    mat	YEn;
    /// observation vector (m*1)
    mat	Yo;
    /// observation ensemble (m*N)
    mat	YoEn;

    /// observation error covariance (m*m)
    mat	R;
    /// error covariance for model forecast (n*n)
    mat	Pf;
    /// error covariance for model analysis (n*n)
    mat	Pa;

    /// Kalman gain (n*m)
    mat	K;

    /// low range of model state (n*1)
    mat	Xl;
    /// high range of model state (n*1)
    mat	Xh;
    /// low range of observation vector (m*1)
    mat	Yl;
    /// high range of observation vector (m*1)
    mat	Yh;

    //unsigned long m_seed;
    bool m_perturb_observe;
    bool m_perturb_model;
};
typedef EnsembleKalmanFilter EnKF;
}
#endif // __LDAS_ENSEMBLEKALMANFILTER_H
