#include "filter.h"

using namespace ldas;
Filter::Filter():state_num(0),observe_num(0)
{
	//ctor
}

Filter::~Filter()
{
	//dtor
}

Filter::Filter(const Filter& other)
{
	//copy ctor
}

Filter::Filter(const unsigned int state,const unsigned int obs):state_num(state),observe_num(obs)
{
	//copy ctor
}
Filter& Filter::operator=(const Filter& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	//assignment operator
	return *this;
}

void Filter::stateNum(const unsigned int n)
{
	state_num=n;
}

unsigned int Filter::stateNum() const
{
	return state_num;
}

void Filter::observeNum(const unsigned int n)
{
	observe_num=n;
}

unsigned int Filter::observeNum() const
{
	return observe_num;
}
