#ifndef __LDAS_QHMODEL_H
#define __LDAS_QHMODEL_H

#include "radiancetransfermodel.h"
#include "permittivity.h"
#include "constant.h"

using namespace std;
namespace ldas
{
///Q-H经验模型，土壤辐射传输模型
///是地表粗糙度的参数化模型,该模型为经验模型,模型采用两个经验性参数 Q 和 h 来描述地表粗糙度对反射率的影响
class QHModel : public RadianceTransferModel
{
public:
    /** Default constructor */
    QHModel():m_q(0),m_h(0) {}
    /** Default destructor */
    virtual ~QHModel() {}
    /** Copy constructor
     *  \param other Object to copy from
     */
    QHModel(const QHModel& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    QHModel& operator=(const QHModel& other);

    QHModel(const SoilParameter& sp);

    virtual void run();
    void soilparameter(const SoilParameter& sp)
    {
        soil=sp;
        m_q=Q(soil.rms,soil.frequency);
        m_h=h(soil.rms,soil.frequency,soil.theta*PI_NUMBER/180);
    }
    SoilParameter soilparameter() const
    {
        return soil;
    }

    Brightness tb() const
    {
        return Tb;
    }
    Emissivity	 emissivity()
    {
        return E;
    }
    Reflectivity GetR()
    {
        return R;
    }
    Emissivity   GetEf()
    {
        return Ef;
    }
    Reflectivity GetRf()
    {
        return Rf;
    }

    /** Get Q value
    * \param rms height of surface roughness (cm) /correlation length of surface (cm)?
    * \param frequency unit Ghz
    * \return Q value
    */
    double Q(const double rms,const double frequency) const;
    double Q() const
    {
        return m_q;
    }
    void Q(const double cq)
    {
    	m_q=cq;
    }
    /** Get h value
    * \param rms height of surface roughness (cm)
    * \param frequency unit Ghz
    * \param theta unit radian
    * \return h value
    */
    double h(const double rms,const double frequency,const double theta) const;
    double h() const
    {
        return m_h;
    }
    void h(const double ch)
    {
    	m_h=ch;
    }
protected:
    SoilParameter soil;
    double m_q,m_h;
private:
    // output
    Brightness		Tb;	// brightness temperature (K)
    Emissivity		E;	// emissivity of rough surface
    Reflectivity	R;	// reflectivity of rough surface
    Emissivity		Ef;	// Fresnel emissivity
    Reflectivity	Rf;	// Fresnel reflectivity

    double GetTao();			// calculate canopy opacity
};

}
#endif // __LDAS_QHMODEL_H
