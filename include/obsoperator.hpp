// obsoperator.hpp
// a general observation operator
// <Xin Li @ CAREERI/CAS> email: lixin@lzb.ac.cn
// August 8, 2005
// Ver 1.0

#ifndef __LDAS_OBSOPERATOR_HPP
#define __LDAS_OBSOPERATOR_HPP

#include "matrix.hpp"

using namespace std;

namespace ldas
{
///观测算子
class ObsOperator
{
public:
    ObsOperator() {}
    virtual ~ObsOperator() {}

    virtual Matrix transfer(const Matrix& X) = 0;
};
}

#endif
