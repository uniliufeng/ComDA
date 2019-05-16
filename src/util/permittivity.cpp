#include "permittivity.h"

using namespace ldas;
Permittivity::Permittivity(const double frequency):m_frequency(frequency)
{
}

Permittivity::~Permittivity()
{
}

double Permittivity::frequency() const
{
	return m_frequency;
}
void Permittivity::frequency(const double frequency)
{
	m_frequency=frequency;
}
// calculate the dielectric constant of free water
complex<double> Permittivity::freewater(const double m_T) const
{
    // parameters and constants
    double Ew0;							// static dielectric constant of water, Ew ~= 80.1
    double Ew00;						// high-frequency limit of Ew0, Ew00 ~= 4.9

    double PI2Tw;						// 2*PI*Tw, Tw: relaxation time of water (s)
    double T	= m_T-273.16;			// temperature of water (Celsius degree)
    double f	= m_frequency*1.0e9;	// electromagnetic frequency (Hz)

    PI2Tw	= 1.1109e-10 - 3.824e-12 * T + 6.938e-14 * T * T - 5.096e-16 * T * T * T;
    Ew0		= 88.045 - 0.4147 * T + 6.295e-4 * T * T + 1.075e-5 * T * T * T;
    Ew00	= 5.27137 + 0.0216474 * T - 0.00131198 * T * T;		// reference ?
    // Ew00	= 4.9;
    complex<double> tmp1(1.0, PI2Tw*f);

    return (Ew00 + (Ew0-Ew00)/tmp1);
}

// calculate the dielectric constant of free water
// reference:	(1) Ulaby et al., 1986. pp2022-2025
complex<double> Permittivity::salinewater(const double m_T, const double m_Salinity) const
{
    double T	= m_T-273.16;			// temperature of water (Celsius degree)
    double f	= m_frequency*1.0e9;	// electromagnetic frequency (Hz)
    double S	= m_Salinity;			// Salinity (¡ë)

    const double E0 = 8.854e-12;		// permitivity of free space (F m^-1)

    double sigma = S*(0.18252 - 1.4619e-3*S + 2.093e-5*S*S - 1.282e-7*S*S*S);
    double delta = 25.0-T;
    double phi = delta*(2.033e-2 + 1.266e-4*delta + 2.464e-6*delta*delta
                        - S*(1.849e-5 - 2.551e-7*delta + 2.551e-8*delta*delta));
    sigma = sigma*exp(-phi);

    complex<double> Efw = freewater(m_T);
    complex<double> tmp = complex<double>(0.0, (sigma/(2*PI*E0*f)) );

    return (Efw-tmp);
}

// calculate the dielectric constant of ice
// reference:	(1) Ulaby et al., 1986. pp2022-2028
//				(2) Kendra et al., 1998
//				(3) Nyfors, 1982
complex<double> Permittivity::ice(const double m_T) const
{
    double T	= m_T;					// temperature of ice (Kelvin)
    double f	= m_frequency*1.0e9;	// electromagnetic frequency (Hz)

    const double real_Eice	= 3.15;
    double image_Eice		= 57.34 * ( 1.0/f + 2.48e-14*sqrt(f) ) * exp(0.0362*T);

    return ( complex<double>(real_Eice, -image_Eice) );
}

// calculate the dielectric constant of dry snow
// reference:	(1) Ulaby et al., 1986. pp2059-2081
//				(2) Matzler, 1996
complex<double> Permittivity::drysnow(const double m_snowDensity,
                                      const double m_T) const
{
    complex<double> Eice = ice(m_T);

    double rs	= m_snowDensity;		// snow density (g cm^-3)

    double real_Edrysnow, image_Edrysnow;

    real_Edrysnow	= 1.0 + 1.5995 * rs + 1.861 * rs * rs * rs;
    image_Edrysnow	= 3.0 * ( rs/0.917 ) * fabs( Eice.imag() )
                      * real_Edrysnow * real_Edrysnow * ( 2.0*real_Edrysnow+1.0 )
                      / ( Eice.real()+2.0*real_Edrysnow )
                      / ( Eice.real()+2.0*real_Edrysnow*real_Edrysnow );

    return ( complex<double>(real_Edrysnow, -image_Edrysnow) );
}

// calculate the dielectric constant of wet snow
// reference:	(1) Ulaby et al., 1986. pp2071-2072

complex<double> Permittivity::wetsnow(const double m_snowDensity,
                                      const double m_snowWetness) const
{
    double rs	= m_snowDensity;		// snow density (g cm^-3)
    double mv	= m_snowWetness;		// snow wetness (%)
    double f	= m_frequency*1.0e9;	// electromagnetic frequency (Hz)

    const double f0 = 9.07e9;
    double A		= 1.0 + 1.83*rs + 0.02*pow(mv, 1.015);
    const double B	= 0.073;
    const double C	= 0.073;
    const double x	= 1.31;
    double tmp		= f/f0;

    double real_Ewetsnow, image_Ewetsnow;

    real_Ewetsnow	= A + B   * pow(mv,x) / (1.0+tmp*tmp);
    image_Ewetsnow	= C * tmp * pow(mv,x) / (1.0+tmp*tmp);

    return ( complex<double>(real_Ewetsnow, -image_Ewetsnow) );
}


// calculate the dielectric constant of soil
// using Dobson mixture model
// Reference:	(1) Ulaby et al., 1986. pp2102
//				(2) Dobson et al., 1985;
complex<double> Permittivity::soil(const double m_LiquidWater,
                                   const double m_T,
                                   const double m_SpecificDensity,
                                   const double m_BulkDensity,
                                   const double m_sand,
                                   const double m_clay) const
{
    complex<double> fw	= freewater(m_T);
    complex<double> ess(4.7, 0.0);		// dielectric constant of solid soil ~= (4.7, 0.0)

    double mv	= m_LiquidWater;		// volumetric water content (m^3 m^-3)
    double ps	= m_SpecificDensity;	// specific density of soil (g cm^-3)
    double pb	= m_BulkDensity;		// bulk density of soil (g cm^-3)
    double S	= m_sand;				// sand fraction (%)
    double C	= m_clay;				// clay fraction (%)

    double a	= 0.65;					// shape factor
    double b	= 1.09-0.11*S/100+0.18*C/100;

    //complex<double> b( (127.48-0.519*S-0.152*C)/100, (133.797-0.603*S-0.166*C)/100 );
    return pow(1.0+(pb/ps)*(pow(ess, a)-1.0)+pow(mv, b)*pow(fw, a)-mv, 1/a);
}

complex<double> Permittivity::soil(const double m_LiquidWater,
                                   const double m_T,
                                   const double m_porosity,
                                   const double m_sand,
                                   const double m_clay) const
{
    complex<double> fw	= freewater(m_T);
    complex<double> ess(4.7, 0.0);		// dielectric constant of solid soil ~= (4.7, 0.0)

    double mv	= m_LiquidWater;		// volumetric water content (m^3 m^-3)
    double n	= m_porosity;			// porosity of soil (-)
    double S	= m_sand;				// sand fraction (%)
    double C	= m_clay;				// clay fraction (%)

    double a	= 0.65;					// shape factor
    double b	= 1.09-0.11*S/100+0.18*C/100;

    return pow(1.0+(1-n)*(pow(ess, a)-1.0)+pow(mv, b)*pow(fw, a)-mv, 1/a);
}

// calculate the penetration depth of a certain media
double Permittivity::depth(complex<double>& permittivity) const
{
    return (LIGHT_SPEED/m_frequency/1.0e9 * sqrt( permittivity.real() )
            / ( 2*PI_NUMBER*fabs(permittivity.imag()) ) );
}
