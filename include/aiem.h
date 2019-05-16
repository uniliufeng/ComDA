#ifndef __LDAS_AIEM_H
#define __LDAS_AIEM_H

#include <iostream>
#include <complex>
#include "radiancetransfermodel.h"
#include "constant.h"

namespace ldas
{
extern "C"
{
    void quagen_(double[], double[], int*);
    void sigma_(double*, double*, double*, double*, double*, double[],
                int*, double*, double*, int*);
}
class AIEM : public RadianceTransferModel
{
public:
    /** Default constructor */
    AIEM():m_mode(AIEM_PASSIVE),m_correlation(AIEM_GAUSSIAN),fresnel_type(AIEM_INCIDENT_ANGLE),npp(128),nss(128) {}
    /** Default destructor */
    virtual ~AIEM() {}
    /** Copy constructor
     *  \param other Object to copy from
     */
    AIEM(const AIEM& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    AIEM& operator=(const AIEM& other);

    AIEM(const SoilParameter& sp);

    // set the corelation type
    void corelation(const CORRELATION_TYPE c=AIEM_GAUSSIAN)
    {
        m_correlation=c;
    }
    CORRELATION_TYPE corelation() const
    {
        return m_correlation;
    }
    // set the approximation type for Fresnel reflection coeff.
    void fresnel(const FRESNEL_TYPE f=AIEM_INCIDENT_ANGLE)
    {
        fresnel_type=f;
    }
    FRESNEL_TYPE fresnel() const
    {
        return fresnel_type;
    }
    // set Gauss integral dimensions
    void gauss_integral(const int m_npp=128, const int m_nss=128);
    virtual void run();
    void active();
    void passive();
    void soilparameter(const SoilParameter& sp)
    {
        soil=sp;
    }
    SoilParameter soilparameter() const
    {
        return soil;
    }
    void permittivity(const std::complex<double> p)
    {
        er=p;
    }
    std::complex<double> permittivity() const
    {
        return er;
    }
    void mode(const AIEM_MODE& m)
    {
        m_mode=m;
    }
    AIEM_MODE mode() const
    {
        return m_mode;
    }
    Emissivity emissivity() const
    {
        return E;    // get the emissivity
    }
    /// 获取后向散射系数，其中h代表hh，v代表vv
    Backscattering scattering() const
    {
        return E;
    }
    Reflectivity GetR()
    {
        return R;    // get the rough surface reflectivity
    }
    Emissivity GetEf()
    {
        return Ef;    // get the Fresnel emissivity
    }
    Reflectivity GetRf()
    {
        return Rf;    // get the Fresnel reflectivity
    }
    Reflectivity GetRco()
    {
        return Rco;    // get the coherent component of R
    }
    Reflectivity GetRnc()
    {
        return Rnc;    // get the non-coherent component of R
    }
protected:
    SoilParameter soil;
    CORRELATION_TYPE	m_correlation;
    FRESNEL_TYPE	fresnel_type;
    std::complex<double>	er;					// permittivity
    AIEM_MODE m_mode;
private:
    // input
    //double			f;					// frequency (GHz)
    //complex<double>	ur;					// permeability
    //double			rms;				// rms height of surface roughness (cm)
    //double			cl;					// correlation length of surface (cm)
    // parameters
    int				npp;				// Gauss integral dimension
    int				nss;				// Gauss integral dimension
    // output
    Emissivity		E;					// emissivity of rough surface
    Reflectivity	R;					// reflectivity of rough surface
    Emissivity		Ef;					// Fresnel emissivity
    Reflectivity	Rf;					// Fresnel reflectivity
    Reflectivity	Rco;				// coherent component of R
    Reflectivity	Rnc;				// non-coherent component of R
};
}
#endif // __LDAS_AIEM_H
