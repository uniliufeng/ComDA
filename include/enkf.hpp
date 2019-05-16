#ifndef __LDAS_ENKF_HPP
#define __LDAS_ENKF_HPP

#include <iostream>
#include <cmath>

#include "matrix.hpp"
#include "random.hpp"
#include "exception.h"

using namespace std;

namespace ldas
{
///Ensemble Kalman filter, 集合卡尔曼滤波算法
class EnKF
{
public:
    ///默认构造函数
    EnKF();

    /** \brief 构造函数
    * \param StateNum 状态数
    * \param ObsNum 观测数
    * \param EnNum 集合数
    */
    EnKF(const unsigned int StateNum, const unsigned int ObsNum, const unsigned int EnNum);
    ///析构函数
    virtual ~EnKF();

    ///参数设置
    void stateNum(const unsigned int num);
    void observeNum(const unsigned int num);
    void ensembleNum(const unsigned int num);
    int stateNum() const;
    int observeNum() const;
    int ensembleNum() const;
    ///method for setting random seed
    void randomSeed(const unsigned long seed);
    ///return random seed
    unsigned long randomSeed() const;

    ///设定模型状态范围
    void xrange(const Matrix<double>& m_Xl, const Matrix<double>& m_Xh);
    ///设定观测范围
    void yrange(const Matrix<double>& m_Yl, const Matrix<double>& m_Yh);
    ///generate model error coveriance matrix Q
    Matrix<double> modelError(const Matrix<double>& m_Xl,const Matrix<double>& m_Xh, const double fac=0.01);
    void is_perturb_model(const bool status);
    bool is_perturb_model() const;

    Matrix<double> observeError(const Matrix<double>& m_Yl,const Matrix<double>& m_Yh, const double fac=0.01);

    ///设置观测
    void observe(const Matrix<double>& m_Yo);
    void is_perturb_observe(const bool status);
    bool is_perturb_observe() const;

    ///预报, generate XfEn, perturb XfEn.
    void predict(const Matrix<double>& m_Q);
    /** update with observation.
    * \param XfEn n*N matrix, model forecast state ensemble
    * \param YEn m*N matrix, ensemble of estimated observation, with H(observe model)
    * \param Yo m*1 matrix, observation
    * \param R m*m, observe error coveriance
    * \param Q n*n, model error coveriance
    */
    void  update(Matrix <double>& XfEn,
                 Matrix <double>& YEn,
                 Matrix <double>& Yo,
                 Matrix <double>& R,
                 Matrix <double>& Q);
    /// update without observation
    /// 返回均值
    void update(Matrix<double>& XfEn);

    ///mean of model state forecast (n*1)
    Matrix<double> GetXf() const;
    ///mean of model state analysis (n*1)
    Matrix<double> GetXa() const;
    ///Matrix holding the Ensemble members forecast (n*N)
    Matrix<double> GetXfEn() const;
    ///Matrix holding the Ensemble members analysis (n*N)
    Matrix<double> GetXaEn() const;
    ///mean of estimated observation (m*1)
    Matrix<double> GetY() const;
    ///ensemble of estimated observation (m*N)
    Matrix<double> GetYEn() const;
    ///observation vector (m*1)
    Matrix<double> GetYo() const;
    ///observation ensemble (m*N)
    Matrix<double> GetYoEn() const;
    ///observation error covariance (m*m)
    Matrix<double> GetR() const;
    ///error covariance for model forecast (n*n)
    Matrix<double> GetPf() const;
    ///error covariance for model analysis (n*n)
    Matrix<double> GetPa() const;
    ///取得增益矩阵K
    Matrix<double> gain() const;
protected:
    ///perturb Yo
    void perturbYo(const Matrix<double>& m_R);

    /// calculate the ensemble mean, retun a matrix
    Matrix<double> ensembleMean(const Matrix<double>& ensemble);
    /// calculate the ensemble covariance
    Matrix<double> ensembleCov(const Matrix<double>& ensemble);
    /// calculate the ensemble perturbation/increment
    Matrix<double> ensembleInc(const Matrix<double>& ensemble);

private:
    /// n: dimension of model state
    unsigned int	n;
    /// m: number of observations
    unsigned int	m;
    /// N: number of emsemble members
    unsigned int	N;

    /// mean of model state forecast (n*1)
    Matrix<double>	Xf;
    /// mean of model state analysis (n*1)
    Matrix<double>	Xa;
    /// Matrix holding the Ensemble members forecast (n*N)
    Matrix<double>	XfEn;
    /// Matrix holding the Ensemble members analysis (n*N)
    Matrix<double>  XaEn;
    /// mean of estimated observation (m*1)
    Matrix<double>	Y;
    /// ensemble of estimated observation (m*N)
    Matrix<double>	YEn;
    /// observation vector (m*1)
    Matrix<double>	Yo;
    /// observation ensemble (m*N)
    Matrix<double>	YoEn;

    /// observation error covariance (m*m)
    Matrix<double>	R;
    /// error covariance for model forecast (n*n)
    Matrix<double>	Pf;
    /// error covariance for model analysis (n*n)
    Matrix<double>	Pa;

    /// Kalman gain (n*m)
    Matrix<double>	K;

    /// low range of model state (n*1)
    Matrix<double>	Xl;
    /// high range of model state (n*1)
    Matrix<double>	Xh;
    /// low range of observation vector (m*1)
    Matrix<double>	Yl;
    /// high range of observation vector (m*1)
    Matrix<double>	Yh;

    unsigned long m_seed;
    bool m_perturb_observe;
    bool m_perturb_model;
};
}

#endif
