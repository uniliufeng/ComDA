#ifndef __LDAS_MEMLS_H
#define __LDAS_MEMLS_H

#include <iostream>
#include <cmath>
#include "radiancetransfermodel.h"
#include "matrix.hpp"
#include "complexmatrix.hpp"
#include "constant.h"

namespace ldas
{
/** \brief MEMLS积雪辐射传输模型.
 The program is about a microwave emission model of layered snowpacks
(MEMLS). The original code is metlab programs, was translated to
C++ by CHE Tao in summer 2003.
When you provide the snowpacks parameters, the program can give you
the brightness temperature in horizatal and vertical polarization
,respectively as the outputs.
The original model was developed by Wiesmann in 1998.
*/
class MEMLS : public RadianceTransferModel
{
public:
    /** Default constructor */
    MEMLS();
    MEMLS(const SnowParameter& sp);
    /** Default destructor */
    virtual ~MEMLS();
    /** Copy constructor
     *  \param other Object to copy from
     */
    MEMLS(const MEMLS& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    MEMLS& operator=(const MEMLS& other);
    virtual void run();
    void parameter(const SnowParameter& sp);
    SnowParameter parameter() const
    {
        return m_snow;
    }
    Brightness brightness() const
    {
        return m_tb;
    }
protected:
    ///calculates the dielectric permittivity of ice
    int epsice();
    ///calculates the dielectric permittivity for dry snow from density.
    void epsr();
    ///fresnel reflection coefficients (assuming eps'' = 0)
    ///(layer n+1 is the air above the snowpack)
    void fresnelrc();
    ///calculates the dielectric permittivity from density for dry snow.
    void ro2epsd();
    ///calculates the permittivity for Wetness > 0
    ///Physical Mixing Model Weise 97 after Matzler 1987 (corrected)
    ///water temperature is assumed constant at 273.15 K
    void mixmod();
    ///computes the absorption coefficient from the dielectric properties
    void abscoeff();
    ///calculates the effective path length in a layer
    void pfadi();
    ///calculates the effective path length in a layer
    ///\param theta incidence angle at snow air interface,unit:radiance
    void pfadc(const double& theta);
    ///calculates the polarization mixing of the interface reflectivities
    ///of each layer (taking into account the first order scattering)
    void polmix();
    ///fresnel reflection coefficients (assuming eps'' = 0)
    ///(layer n+1 is the air above the snowpack)
    void fresnelc();
    ///locates and treats coherent layers in a snowpack
    void slred();
    ///calculates the scattering coefficient from structural parameters
    ///different algorithms can be chosen, by changing "sccho"
    ///Note: this sccho=11 is used as snow layers
    void sccoeff(
        Matrix <double>& gbih,
        Matrix <double>& gbiv,
        Matrix <double>& ga2i);
    ///calculates the layer reflectivity and transmissivity
    void rt(
        Matrix <double>& gbi		//gbi:  scattering coefficient
    );
    ///calculates the upwelling brightness temperatures D
    ///si:   interface reflectivity
    void layer(Matrix <double>& si);
private:
    ///eice:  dielectric permittivity of ice
    Matrix <double> eice;
    Matrix <double> gb6;
    Matrix <double> gc6;
    Matrix <double> gf6;
    ///gs6:  6-flux scattering coefficient
    Matrix <double> gs6;

    ///epsi: real part of dielectric permittivity
    Matrix <double> epsi;
    ///epsii: imaginary part of dielectric permittivity
    Matrix <double> epsii;

    ///dei:  effective path length [m]
    Matrix <double> dei;
    ///tei:  local incidence angle
    Matrix <double> tei;
    ///tscat: tau scat
    Matrix <double> tscat;

    ///sih:  interface reflectivity at h pol
    Matrix <double> sih;
    ///siv:  interface reflectivity at v pol
    Matrix <double> siv;

    ///gbih:  2-flux scattering coefficient at h pol
    Matrix <double> gbih;
    ///gbiv:  2-flux scattering coefficient at v pol
    Matrix <double> gbiv;
    ///ga2i:  2-flux absorption coefficient
    Matrix <double> ga2i;

    ///gai:  absorption coefficient
    Matrix <double> gai;

    ///ri:   layer reflectivity
    Matrix <double> ri;
    ///ti:   layer transmissivity
    Matrix <double> ti;

    ///FH:   Fresnel reflection coefficient at h pol
    Matrix <double> FH;
    ///FV:   Fresnel reflection coefficient at v pol
    Matrix <double> FV;

    ///D:    upwelling brightness temperature
    Matrix <double> D;
    ///------------replace with m_snow class
    //int col;//雪的层数
    ///雪层的编号
    Matrix <int> num;
    //Matrix <double> Ti;//每层的温度         (K)
    //Matrix <double> Wi;///每层的湿度        (g/cm^2)
    //Matrix <double> roi;//每层的密度        (g/cm^3)
    //Matrix <double> di;//每层的厚度         (cm)
    //Matrix <double> pci;//每层的相关长度     (m)
    //以上为memls类的参数
protected:
    ///模型计算结果：亮度温度
    Brightness m_tb;
    ///模型参数
    SnowParameter m_snow;
};
}
#endif // __LDAS_MEMLS_H
