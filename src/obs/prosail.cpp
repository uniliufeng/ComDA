#include "prosail.h"

using namespace ldas;
ProSail::ProSail()
{
	//ctor
}

ProSail::~ProSail()
{
	//dtor

}
ProSail::ProSail(const CanopyParameter& cp)
{
	m_canopy=cp;
}
ProSail::ProSail(const ProSail& other)
{
	//copy ctor
	m_canopy=other.parameter();
}

ProSail& ProSail::operator=(const ProSail& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	//assignment operator
	m_canopy=rhs.parameter();
	return *this;
}
void ProSail::run()
{
	pro4sail_(&m_canopy.chlorophyll_content,&m_canopy.carotenoid_content,&m_canopy.brownpigment_content,
			&m_canopy.cw,&m_canopy.cm,&m_canopy.structure_coefficient, &m_canopy.angle, &m_canopy.lai,&m_canopy.hot_spot,
			&m_canopy.solar_zenith,&m_canopy.observer_zenith,&m_canopy.azimuth,&m_canopy.soil_coefficient,
			&m_canopy.skyl,resh,resv);
}
