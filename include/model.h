#ifndef __LDAS_BaseModel_H
#define __LDAS_BaseModel_H

#include "time.hpp"
namespace ldas
{
class BaseModel
{
public:
    /** Default constructor */
    BaseModel();
    /** Default destructor */
    virtual ~BaseModel();
    /** Copy constructor
     *  \param other Object to copy from
     */
    BaseModel(const BaseModel& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    BaseModel& operator=(const BaseModel& other);

    /** Virtual method, full parameter configuration */
    virtual void run()=0;
    /** Run model with single point(index)
    * \param spatial_index 1D array position
    */
    virtual void run(const int spatial_index) {}
    virtual void run(const Time& time_point) {}

    Time time() const
    {
        return m_time;
    }
    void time(const Time& t)
    {
        m_time=t;
    }
protected:
    ///current BaseModel time
    Time m_time;
};

class ModelConfig
{
public:
    ModelConfig()
    {
        time_step=3600; // 1 hour
        output_span=1;
        restart_span=240;//10 days
        is_regional=true;
        total_patches=0;
        total_steps=0;
    }
    virtual ~ModelConfig() {}
    ///if true, the path is a directory, else it is a filename
    bool is_regional;
    std::string forcing_path,vegetation_path,soil_path;
    std::string restart_path,output_path;
    ///unit: time_step
    int restart_span;
    int output_span;
    ///unit:seconds
    int time_step;
    ///model step for simulation
    int total_steps;
    ///total patches
    int total_patches;
};
class ModelParameter
{
public:
    ModelParameter():numpatch(0) {}
    ModelParameter(const int p):numpatch(p) {}
    virtual ~ModelParameter() {}
    virtual int patches() const
    {
        return numpatch;
    }
    virtual void patches(const int p)
    {
        numpatch=p;
    }
protected:
    int numpatch;
};
class ModelForcing
{
public:
    ModelForcing() {}
    virtual ~ModelForcing() {}

    /** \brief 计算可见和近红波段的直接和散射辐射, according to SiB2 model code
    * \param jday Julian cal day (1..365)
    * \param lat  Centered latitude (radians)
    * \param lon  Centered longitude (radians)
    * \param swdown downward short wave radiation
    * \param sols atm vis direct beam solar rad onto srf [W/m2] (output variable)
    * \param soll atm nir direct beam solar rad onto srf [W/m2] (output variable)
    * \param solsd atm vis diffuse solar rad onto srf [W/m2] (output variable)
    * \param solld atm nir diffuse solar rad onto srf [W/m2] (output variable)
    * \return void
    */
    void radiation(int jday,int msec,double lon,double lat, double swdown, double &sols, double &soll, double &solsd, double &solld);

    /// longwave radiation downward (W m^-2)
    double lwdown;
    /// shortwave radiation downward (W m^-2)
    double swdown;
    /// air temperature (K)
    double airtemperature;
    /// air pressure (mb)
    double airpressure;
    /// specific humidity (kg kg^-1)
    double humidity;
    /// rainfall? (mm)?
    double precipitation;
    /// wind speed U (m s^-1)
    double uwind;
    /// wind speed V (m s^-1)
    double vwind;
};
}
#endif // __LDAS_BaseModel_H
