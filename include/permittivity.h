#ifndef __LDAS_PERMITTIVITY_H
#define __LDAS_PERMITTIVITY_H
#include "constant.h"
#include <cmath>
#include <complex>

using namespace std;

namespace ldas
{
class Permittivity
{
public:
    /** Default constructor
    	\param frequency 频率，单位：GHz
    */
    Permittivity(const double frequency);
    /** Default destructor */
    virtual ~Permittivity();

    double frequency() const;
    void frequency(const double frequency);

    /** \brief calculate the dielectric constant of free water

    	reference 1:
    	Ulaby, F.T., Moore, R.K., & Fung, A.K. 1986. Microwave Remote Sensing:
    	Active and Passive. Volume III, from Theory to Applications.
    	Artech House, Norwood, MA.

    	reference 2:
    	Dobson, M.C., Ulaby, F.T., Hallikainen, M.T., & El-Rayes, M.A., 1985.
    	Microwave dielectric behavior of wet soil: part II: dielectric mixing models.
    	IEEE Transactions on Geoscience and Remote Sensing GE-23 (1), 35-46.
    	\param m_T 温度，单位：K
    */
    complex<double> freewater(const double m_T) const;
    /** \brief calculate the dielectric constant of free water

     reference:	(1) Ulaby et al., 1986. pp2022-2025
    	\param m_T 温度，单位：K
    	\param m_frequency 频率，单位：GHz
    	\param m_Salinity 盐度，盐分，盐浓度，单位：？
    */
    complex<double> salinewater(const double m_T,
                                const double m_Salinity) const;

    /** \brief calculate the dielectric constant of ice

        	reference 1:
        	Ulaby, F.T., Moore, R.K., & Fung, A.K. 1986. Microwave Remote Sensing:
        	Active and Passive. Volume III, from Theory to Applications.
        	Artech House, Norwood, MA.

    	reference 2:
    	Kendra, J.R., Sarabandi, K., & Ulaby, F.T. 1998.
    	Radar measurements of snow: experiment and analysis.
    	IEEE Transactions on Geoscience and Remote Sensing, 36: 864-879

    	reference 3:
    	Nyfors, E., 1982. On the dielectric properties of dry snow in
    	the 800 MHz to 13 GHz range. Radio Lab.,
    	Helsinki Univ. Technol., Helsinki, Finland, Tech. Rep. S135
    	\param m_T 温度，单位：K
    */
    complex<double> ice(const double m_T) const;

    /** \brief calculate the dielectric constant of dry snow

    	reference 2:
    	Matzler, C. 1996. Microwave permittivity of dry snow.
    	IEEE Transactions on Geoscience and Remote Sensing, 34, 573-581
    	\param m_snowDensity 雪密度，单位：?
        \param m_T 温度，单位：K
    */
    complex<double> drysnow(const double m_snowDensity,
                            const double m_T) const;
    /** \brief calculate the dielectric constant of wet snow

     reference:	(1) Ulaby et al., 1986. pp2071-2072
    \param m_snowDensity 雪密度，单位：?
    \param m_snowWetness 雪湿度，单位：
     */
    complex<double> wetsnow(const double m_snowDensity,
                            const double m_snowWetness) const;

    /** \brief calculate the dielectric constant of soil.

     using Dobson mixture model
     Reference:	(1) Ulaby et al., 1986. pp2102
    			(2) Dobson et al., 1985;
    */
    complex<double> soil(const double m_LiquidWater,
                         const double m_T,
                         const double m_SpecificDensity,
                         const double m_BulkDensity,
                         const double m_sand,
                         const double m_clay) const;
    /// use soil porosity as input
    complex<double> soil(const double m_LiquidWater,
                         const double m_T,
                         const double m_porosity,
                         const double m_sand,
                         const double m_clay) const;
    /** calculate the penetration depth of a certain media
    *
    * \note conditions: p.imag/p.real<0.1
    */
    double depth(complex<double>& permittivity) const;
private:
    double m_frequency;
};
}
#endif // __LDAS_PERMITTIVITY_H


