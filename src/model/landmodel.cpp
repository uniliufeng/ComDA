#include "landmodel.h"

using namespace ldas;
LandModel::LandModel()
{
	//ctor
}

LandModel::~LandModel()
{
	//dtor
}

LandModel::LandModel(const LandModel& other)
{
	//copy ctor
}

LandModel& LandModel::operator=(const LandModel& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	//assignment operator
	return *this;
}
