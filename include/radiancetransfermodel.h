#ifndef __LDAS_RADIANCETRANSFERMODEL_H
#define __LDAS_RADIANCETRANSFERMODEL_H
#include <iostream>
#include "matrix.hpp"
namespace ldas
{
///radiance transfer model, the abstract class
class RadianceTransferModel
{
public:
    /** Default constructor */
    RadianceTransferModel() {}
    /** Default destructor */
    virtual ~RadianceTransferModel() {}
    /** Copy constructor
     *  \param other Object to copy from
     */
    RadianceTransferModel(const RadianceTransferModel& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    RadianceTransferModel& operator=(const RadianceTransferModel& other);
    virtual void run()=0;
protected:
private:
};
///define soil parameter required by soil radiance transfer model
class SoilParameter
{
public:
    SoilParameter();
    virtual ~SoilParameter() {}
    SoilParameter& operator=(const SoilParameter& other);
    SoilParameter(const SoilParameter& other);
    bool has_vegetation() const
    {
        return m_vegetated;
    }
    void has_vegetation(const bool vegetated)
    {
        m_vegetated=vegetated;
    }
    // sensor parameters
    double theta;		// viewing angle of sensor (degree)
    double frequency;	// frequency (GHz)
    // vegetation parameters
    double Op;			// single scattering albedo of vegetation (-)
    double b1;			// opacity coefficient (-)
    double x;			//
    // soil parameters
    double rms;		// rms height of surface roughness (cm)
    double cl;		// correlation length of surface (cm)
    //unused in Qh?
    double ps;			// specific density of soil (g cm^-3)
    //unused in Qh
    double pb;			// bulk density of soil (g cm^-3)
    double sand;		// % Sand
    double clay;		// % Clay
    double porosity;
    // media variables of land surface
    double Me;			// surface soil moisture (m-3/m-3)
    double We;			// vegetation water content (kg m^-2)
    double Te;			// surface temperature (K)
private:
    bool m_vegetated;
};
///define snow parameter required by snow radiance transfer model
class SnowParameter
{
public:
    SnowParameter();
    virtual ~SnowParameter() {}
    SnowParameter& operator=(const SnowParameter& other);
    SnowParameter(const SnowParameter& other);

    int    col; //number of snow layers
    //sensor parameter
    double frequency;               //frequnecy (GHz)
    double theta;                   //Nadir angle (Degrees)
    double temperature_sky;                    //brightness temperature of sky (K)
    double temperature_ground;                    //temperature of ground (K)

    ///\todo: convert to dynamic array?

    //Matrix<int>    num;		    	//number of layers
    ///temperature of each layer [K]
    Matrix<double> temperature;
    ///wetness of each layer [g/cm2]
    Matrix<double> wetness;
    ///ice content of each layer [g/cm2]
    Matrix<double> icecontent;
    ///density of each layer [g/cm3]
    Matrix<double> density;
    ///depth of each layer   [cm]
    Matrix<double> depth;
    ///correction length     [mm]
    Matrix<double> pci;
};
class CanopyParameter
{
public:
    CanopyParameter();
    virtual ~CanopyParameter() {}
    CanopyParameter& operator=(const CanopyParameter& other);
    CanopyParameter(const CanopyParameter& other);
    ///\brief angle parameters
    ///solar zenith angle (°) 太阳天顶角
    double solar_zenith;
    ///observer zenith angle (°) 视场天顶角
    double observer_zenith;
    ///azimuth (°) 方位角
    double azimuth;

    ///leaf parameters
    ///chlorophyll content (ug.cm-2) 叶绿素含量
    double chlorophyll_content;
    ///carotenoid content (ug.cm-2) 类胡萝卜素含量
    double carotenoid_content;
    ///brown pigment content (arbitrary units) 褐色色素含量
    double brownpigment_content;
    ///EWT (cm) 等效水厚度
    double cw;
    ///LMA (g.cm-2)?
    double cm;

    ///canopy parameters
    ///structure coefficient 结构系数
    double structure_coefficient;
    ///leaf area index 叶面积指数
    double lai;
    ///average leaf angle 平均叶角度,unit:degree
    double angle;

    ///soil parameter
    ///soil coefficient 土壤系数
    double soil_coefficient;

    ///% diffuse/direct radiation 漫射/直射辐射
    double skyl;
    ///hot spot 热点
    double hot_spot;
private:
};
class Brightness
{
public:
    double h;
    double v;
    Brightness (const double h=0.0, const double v=0.0):h(h),v(v)
    {
    }
    virtual ~Brightness() {}
    Brightness& operator=(const Brightness& other);
    Brightness(const Brightness& other);
    friend std::ostream & operator<<(std::ostream& os,const Brightness& b)
    {
        os<<"("<<b.h<<","<<b.v<<")";
        return os;
    }
};

typedef Brightness Reflectivity;
typedef Brightness Emissivity;
typedef Brightness Backscattering;
}
#endif // RADIANCETRANSFERMODEL_H
