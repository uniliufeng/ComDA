#include "radiancetransfermodel.h"

using namespace ldas;

RadianceTransferModel::RadianceTransferModel(const RadianceTransferModel& other)
{
//copy ctor
}

RadianceTransferModel& RadianceTransferModel::operator=(const RadianceTransferModel& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

SoilParameter::SoilParameter():sand(0),clay(0),m_vegetated(false),porosity(0),Me(0),We(0),Te(0),ps(0),pb(0),x(0),b1(0),Op(0),theta(0),frequency(0),cl(0),rms(0)
{

}

SoilParameter::SoilParameter(const SoilParameter& other)
{
    m_vegetated=other.has_vegetation();
    sand=other.sand;
    clay=other.clay;
    porosity=other.porosity;
    Me=other.Me;
    We=other.We;
    Te=other.Te;
    ps=other.ps;
    pb=other.pb;
    x=other.x;
    b1=other.b1;
    Op=other.Op;
    theta=other.theta;
    frequency=other.frequency;
    rms=other.rms;
    cl=other.cl;
}

SoilParameter& SoilParameter::operator=(const SoilParameter& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    m_vegetated=other.has_vegetation();
    sand=other.sand;
    clay=other.clay;
    porosity=other.porosity;
    Me=other.Me;
    We=other.We;
    Te=other.Te;
    ps=other.ps;
    pb=other.pb;
    x=other.x;
    b1=other.b1;
    Op=other.Op;
    theta=other.theta;
    frequency=other.frequency;
    rms=other.rms;
    cl=other.cl;
    return *this;
}

SnowParameter::SnowParameter():col(0),frequency(0),theta(0),temperature_sky(0),temperature_ground(0)
{

}
SnowParameter::SnowParameter(const SnowParameter& other)
{
    col=other.col;
    frequency=other.frequency;
    theta=other.theta;
    temperature_sky=other.temperature_sky;
    temperature_ground=other.temperature_ground;
    if (col>0)
    {
        temperature=other.temperature;
        wetness=other.wetness;
        icecontent=other.icecontent;
        density=other.density;
        depth=other.depth;
        pci=other.pci;
    }
}
SnowParameter& SnowParameter::operator=(const SnowParameter& other)
{
    if (this == &other) return *this; // handle self assignment
    col=other.col;
    frequency=other.frequency;
    theta=other.theta;
    temperature_sky=other.temperature_sky;
    temperature_ground=other.temperature_ground;
    if (col>0)
    {
        temperature=other.temperature;
        wetness=other.wetness;
        icecontent=other.icecontent;
        density=other.density;
        depth=other.depth;
        pci=other.pci;
    }
    return *this;
}

CanopyParameter::CanopyParameter():solar_zenith(0),observer_zenith(0),azimuth(0),chlorophyll_content(0),
    carotenoid_content(0),brownpigment_content(0),cw(0),cm(0),structure_coefficient(0),lai(0),angle(0),soil_coefficient(0),
    skyl(0),hot_spot(0)
{

}
CanopyParameter::CanopyParameter(const CanopyParameter& other)
{
	solar_zenith=other.solar_zenith;
	observer_zenith=other.observer_zenith;
	azimuth=other.azimuth;
	chlorophyll_content=other.chlorophyll_content;
    carotenoid_content=other.carotenoid_content;
    brownpigment_content=other.brownpigment_content;
    cw=other.cw;
    cm=other.cm;
    structure_coefficient=other.structure_coefficient;
    lai=other.lai;
    angle=other.angle;
    soil_coefficient=other.soil_coefficient;
    skyl=other.skyl;
    hot_spot=other.hot_spot;
}
CanopyParameter& CanopyParameter::operator=(const CanopyParameter& other)
{
    if (this == &other) return *this; // handle self assignment
	solar_zenith=other.solar_zenith;
	observer_zenith=other.observer_zenith;
	azimuth=other.azimuth;
	chlorophyll_content=other.chlorophyll_content;
    carotenoid_content=other.carotenoid_content;
    brownpigment_content=other.brownpigment_content;
    cw=other.cw;
    cm=other.cm;
    structure_coefficient=other.structure_coefficient;
    lai=other.lai;
    angle=other.angle;
    soil_coefficient=other.soil_coefficient;
    skyl=other.skyl;
    hot_spot=other.hot_spot;
    return *this;
}

Brightness& Brightness::operator=(const Brightness& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    h=other.h;
    v=other.v;
    return *this;
}

Brightness::Brightness(const Brightness& other)
{
    h=other.h;
    v=other.v;
}
