// modeloperator.cpp
// implementation of a general model operator
// <Xin Li @ CAREERI/CAS> email: lixin@lzb.ac.cn
// August 8, 2005
// Dec 5, 2005
// Ver 1.0

#include "modeloperator.hpp"

using namespace ldas;

ModelOperator::ModelOperator()
{
}

ModelOperator::~ModelOperator()
{
}

int ModelOperator::GetCurrentStep() const
{
    return CurrentStep;
}

Time ModelOperator::GetCurrentTime() const
{
    return CurrentTime;
}

void ModelOperator::SetEnsembleNum(const int m_EnsembleNum)
{
    EnsembleNum = m_EnsembleNum;
}

