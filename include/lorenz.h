#ifndef __LDAS_LORENZ_H
#define __LDAS_LORENZ_H

#include "model.h"

namespace ldas
{
///Lorenz Model, 可用于测试数据同化算法
/** The Lorenz equations (lorenz.h) */
/**
* This header file describes function type of the Lorenz equations
* for the Runge-Kutta routine (rk.h).
* The Lorenz equations derived by simplifying of convection rolls arising
* in the atmosphere, and described as ordinary differential equations
* that have 3 variables:
* dx / dt = -δx + δy,
* dy / dt = -xz + rx - y,
* dz / dt = xy - bz,
* , where δ, b, and r are parameters, and this time, they are set to as follows:
* δ = 10, b = 8 / 3, and r = 28.
* The Lorenz equations exhibit sensitive dependence on initial conditions.
* This phenomenon is used as the root of "the butterfly effectexternal link".
*/
class Lorenz : public BaseModel
{
public:
    /** Default constructor */
    Lorenz();
    Lorenz(const double sigma, const double lambda, const double beta);
    /** Default destructor */
    virtual ~Lorenz();
    /** Copy constructor
     *  \param other Object to copy from
     */
    Lorenz(const Lorenz& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    Lorenz& operator=(const Lorenz& other);
    virtual void run();
    void parameter(const double X[3]);
    const double* parameter() const;
    const double* result() const;
    void output() const;
protected:
    double m_sigma,m_lambda,m_beta,m_h;
    double m_x[3],m_x0[3];
    /// Return the derivatives (dx/dt dy/dt dz/dt)
    void lorzrk(const double X[3], double deriv[3]);
    /**
    ** Putting the Runge-Kutta step 'h,' the rank of the differential equation N,
    ** and the vector 'X0' which contains the state variable at time 't'
    ** to get the state variable 'X1' at time 't + h'.
    */
    void rk4(double h, double X0[3], double X[3]); /* Runge-Kutta */
private:

};
}
#endif // __LDAS_LORENZ_H
