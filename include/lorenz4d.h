#ifndef __LDAS_Lorenz4d_H
#define __LDAS_Lorenz4d_H

#include "model.h"
#include <armadillo>

using namespace arma;
namespace ldas
{
///Lorenz4d Model, 可用于测试数据同化算法
/** The Lorenz4d equations (Lorenz4d.h) */
/**
* This header file describes function type of the Lorenz4d equations
* for the Runge-Kutta routine (rk.h).
* The Lorenz4d equations derived by simplifying of convection rolls arising
* in the atmosphere, and described as ordinary differential equations
* that have 3 variables:
* dx / dt = -δx + δy,
* dy / dt = -xz + rx - y,
* dz / dt = xy - bz,
* , where δ, b, and r are parameters, and this time, they are set to as follows:
* δ = 10, b = 8 / 3, and r = 28.
* The Lorenz4d equations exhibit sensitive dependence on initial conditions.
* This phenomenon is used as the root of "the butterfly effectexternal link".
*/
class Lorenz4d : public BaseModel
{
public:
    /** Default constructor */
    Lorenz4d();
    Lorenz4d(const double F, const double h, const unsigned int n=40);
    /** Default destructor */
    virtual ~Lorenz4d();
    /** Copy constructor
     *  \param other Object to copy from
     */
    Lorenz4d(const Lorenz4d& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Lorenz4d& operator=(const Lorenz4d& other);
    virtual void run();
    void parameter(const mat& X);
    mat parameter() const;
    mat result() const;
    void output() const;
protected:
    double m_h,m_F;
    unsigned int m_nvar; ///X[j]
    mat m_X;
    /// Return the derivatives (dx/dt dy/dt dz/dt)
    void lorzrk(const mat& X, const double F, mat& deriv);
    /**
    ** Putting the Runge-Kutta step 'h,' the rank of the differential equation N,
    ** and the vector 'X0' which contains the state variable at time 't'
    ** to get the state variable 'X1' at time 't + h'.
    */
    void rk4(const double h,const double F); /* Runge-Kutta */
private:

};
}
#endif // __LDAS_Lorenz4d_H
