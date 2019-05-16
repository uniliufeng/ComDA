// modeloperator.hpp
// a general model operator
// <Xin Li @ CAREERI/CAS> email: lixin@lzb.ac.cn
// August 8, 2005
// October 13, 2005
// Dec 5, 2005
// Ver 1.0

#ifndef __LDAS_MODELOPERATOR_HPP
#define __LDAS_MODELOPERATOR_HPP

#include "matrix.hpp"
#include "time.hpp"

namespace ldas
{
///模型算子
class ModelOperator
{
protected:
    /// modeling begin time
    Time	BeginTime;
    /// modeling end time
    Time	EndTime;
    /// current time
    Time	CurrentTime;
    /// model step (seconds)
    double	DeltaModel;
    /// ensemble number
    unsigned int EnsembleNum;
    /// total steps of model
    unsigned int TotalSteps;
    /// current step in modeling
    unsigned int CurrentStep;
public:
    ///默认构造函数
    ModelOperator();
    ///析构函数
    virtual ~ModelOperator();
    ///获取当前运行到的模型位置
    int GetCurrentStep() const;
    ///获取当前模型计算到的时间点
    Time GetCurrentTime() const;
    ///设置集合数
    void SetEnsembleNum(const int m_EnsembleNum);

    virtual void Init() = 0;
    virtual void Predict() = 0;
    ///计算模型下一步或下几步
    //virtual Matrix step(const Matrix& Xa, const int steps) = 0;
};
}

#endif
