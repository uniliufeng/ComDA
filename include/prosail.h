#ifndef __LDAS_PROSAIL_H
#define __LDAS_PROSAIL_H

#include "radiancetransfermodel.h"

namespace ldas
{
extern "C"
{
    void pro4sail_(double*,double*,double*,double*,double*,double*,double*,
                   double*,double*,double*,double*,double*,double*,double*,double[],double[]);
}
class ProSail : public RadianceTransferModel
{
public:
    /** Default constructor */
    ProSail();
    ProSail(const CanopyParameter& cp);
    /** Default destructor */
    virtual ~ProSail();
    /** Copy constructor
     *  \param other Object to copy from
     */
    ProSail(const ProSail& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    ProSail& operator=(const ProSail& other);
    virtual void run();
    void parameter(const CanopyParameter& cp)
    {
        m_canopy=cp;
    }
    CanopyParameter parameter() const
    {
        return m_canopy;
    }
    void output()
    {
    	for(int i=0;i<2101;i++)
		{
			std::cout<<i+400<<" "<<resh[i]<<" "<<resv[i]<<std::endl;
		}
    }
protected:
    CanopyParameter m_canopy;
    ///半球反射率
    double resh[2101];
    ///方向性反射率
    double resv[2101];
private:
};
}
#endif // __LDAS_PROSAIL_H
