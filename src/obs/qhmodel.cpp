#include "qhmodel.h"

using namespace ldas;

QHModel::QHModel(const QHModel& other)
{
    //copy ctor
    soil=other.soilparameter();
	m_q=Q(soil.rms,soil.frequency);
	m_h=h(soil.rms,soil.frequency,soil.theta*PI_NUMBER/180);
}

QHModel& QHModel::operator=(const QHModel& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    soil=rhs.soilparameter();
	m_q=Q(soil.rms,soil.frequency);
	m_h=h(soil.rms,soil.frequency,soil.theta*PI_NUMBER/180);
    return *this;
}

QHModel::QHModel(const SoilParameter& sp):soil(sp)
{
	m_q=Q(soil.rms,soil.frequency);
	m_h=h(soil.rms,soil.frequency,soil.theta*PI_NUMBER/180);
}
// forward model to calculate the brightness temperature
void QHModel::run()
{
    Permittivity p(soil.frequency);
    complex<double> er	= p.soil(soil.Me, soil.Te, soil.porosity, soil.sand, soil.clay);

    // calculate the Fresnel reflectivity
    double theta=soil.theta*PI_NUMBER/180.0;
    double cs		= std::cos(theta);
    double si2		= std::sin(theta) * std::sin(theta);
    complex<double> stem	= sqrt(er-si2);
    Rf.h	= abs( (cs-stem)/(cs+stem) );
    Rf.v	= abs( (er*cs-stem)/(er*cs+stem) );
    Rf.h	= Rf.h * Rf.h;
    Rf.v	= Rf.v * Rf.v;
    Ef.h	= 1.0-Rf.h;
    Ef.v	= 1.0-Rf.v;

    // calculate the reflectivity of rough soil
    R.h		= ((1-m_q)*Rf.h+m_q*Rf.v) * m_h;
    R.v		= ((1-m_q)*Rf.v+m_q*Rf.h) * m_h;
    E.h		= 1.0-R.h;
    E.v		= 1.0-R.v;

    if (soil.has_vegetation())
    {
        double tr	= std::exp( -GetTao() );		// transmissivity
        Tb.h	= soil.Te * ( E.h*tr+(1-soil.Op)*(1-tr)*(1+R.h*tr) );
        Tb.v	= soil.Te * ( E.v*tr+(1-soil.Op)*(1-tr)*(1+R.v*tr) );
    }
    else
    {
        Tb.h	= soil.Te*E.h;
        Tb.v	= soil.Te*E.v;
    }
}

// calculate the opacity of vegetation
// reference: Jackson and Schmugge, 1991
double QHModel::GetTao()
{
    // "b" depends on canopy structure and frequency
    // the unit of wavelength is transfered to cm
    double b;	// opacity coefficient (-)
    double wavelength = LIGHT_SPEED/soil.frequency/1.0e9;

    b	= soil.b1 * pow(wavelength*100, soil.x);

    return ( b * soil.We/std::cos(soil.theta*PI_NUMBER/180) );
}


// references:
// 1. Wang and Choudhury, 1981
double QHModel::Q(const double rms,const double frequency) const
{
    return 0.35*( 1.0-std::exp(-0.6*rms*rms*frequency) );
}

double QHModel::h(const double rms,const double frequency,const double theta) const
{
    //double theta=soil.theta*PI_NUMBER/180;
    double wavenumber	= 2.0*PI_NUMBER*frequency*1.0e9/LIGHT_SPEED;	// m^-1
    double ks		= (rms*0.01)*wavenumber;					// wavenumber -> cm^-1

    return std::exp( -(2.0*ks*std::cos(theta))*(2.0*ks*std::cos(theta)) );
}

// [1] Wang JR, Choudhury BJ. Remote sensing of soil moisture content over bare
// field at 1.4 GHz frequency[J]. Journal of Geophysical Research-oceans
// and Atmospheres, 1981, 86(NC6): 5277£­5282.
