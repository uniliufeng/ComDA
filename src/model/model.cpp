#include "model.h"

using namespace ldas;

BaseModel::BaseModel()
{
    //ctor
}

BaseModel::~BaseModel()
{
    //dtor
}

BaseModel::BaseModel(const BaseModel& other)
{
    //copy ctor
}

BaseModel& BaseModel::operator=(const BaseModel& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
