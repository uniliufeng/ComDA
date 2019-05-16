#ifndef LANDMODEL_H
#define LANDMODEL_H

#include "model.h"
#include "time.hpp"
namespace ldas
{

class LandModel : public Model
{
public:
    /** Default constructor */
    LandModel();
    /** Default destructor */
    virtual ~LandModel();
    /** Copy constructor
     *  \param other Object to copy from
     */
    LandModel(const LandModel& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    LandModel& operator=(const LandModel& other);
    /** Access m_time
     * \return The current value of m_time
     */
    Time Gettime()
    {
        return m_time;
    }
    /** Set m_time
     * \param val New value to set
     */
    void Settime(Time val)
    {
        m_time = val;
    }
    /** Access m_longitude
     * \return The current value of m_longitude
     */
    double Getlongitude()
    {
        return m_longitude;
    }
    /** Set m_longitude
     * \param val New value to set
     */
    void Setlongitude(double val)
    {
        m_longitude = val;
    }
    /** Access m_lattitude
     * \return The current value of m_lattitude
     */
    double Getlattitude()
    {
        return m_lattitude;
    }
    /** Set m_lattitude
     * \param val New value to set
     */
    void Setlattitude(double val)
    {
        m_lattitude = val;
    }

    virtual void run()=0;
protected:
    Time m_begintime;
    TimeSpan m_timedelta;
    Time m_endtime;
    double m_longitude;
    double m_lattitude;
private:
};

}
#endif // LANDMODEL_H
