#include "lorenz.h"

using namespace ldas;
Lorenz::Lorenz():m_sigma(10.0),m_lambda(28),m_beta(8.0/3),m_h(0.01)
{
    //ctor
}

Lorenz::Lorenz(const double sigma, const double lambda, const double beta)
    :m_sigma(sigma),m_lambda(lambda),m_beta(beta),m_h(0.01)
{

}

Lorenz::~Lorenz()
{
    //dtor
}

Lorenz::Lorenz(const Lorenz& other)
{
    //copy ctor
    m_sigma=10;
    m_lambda=28;
    m_beta=8.0/3;
    m_h=0.01;
    for(int i=0; i<3; i++)
    {
        m_x0[i]=other.parameter()[i];
    }
}

Lorenz& Lorenz::operator=(const Lorenz& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    m_sigma=10;
    m_lambda=28;
    m_beta=8.0/3;
    m_h=0.01;
    for(int i=0; i<3; i++)
    {
        m_x0[i]=other.parameter()[i];
    }
    return *this;
}

void Lorenz::parameter(const double X[3])
{
    for(int i=0; i<3; i++)
        m_x0[i]=X[i];
}
const double* Lorenz::parameter() const
{
    return m_x0;
}
const double* Lorenz::result() const
{
    return m_x;
}

void Lorenz::run()
{
    rk4(m_h,m_x0,m_x);
}

void Lorenz::rk4(double h,double X0[3], double X[3])
{
    double d1[3],d2[3],d3[3],dX[3],Xa[3];

    /* d1 = hF(x(t), t) */
    lorzrk(X0, dX);
    for(int i = 0; i<3; i++)
    {
        d1[i] = h * dX[i];
        Xa[i] = X0[i] + 0.5 * d1[i];
    }

    /* d2 = hF(x(t) + d1 / 2, t + h / 2) */
    lorzrk(Xa, dX);
    for(int i = 0; i < 3; i++)
    {
        d2[i] = h * dX[i];
        Xa[i] = X0[i] + 0.5 * d2[i];
    }

    /* d3 = hF(x(t) + d2 / 2, t + h / 2) */
    lorzrk(Xa, dX);
    for(int i = 0; i < 3; i++)
    {
        d3[i] = h * dX[i];
        Xa[i] = X0[i] + d3[i];
    }

    /* x(t + h) = x(t) + d1 / 6 + d2 / 3 + d3 / 3 */
    lorzrk(Xa, dX);
    for(int i = 0; i < 3; i++)
        X[i] = X0[i] + (d1[i] + d2[i] * 2 + d3[i] * 2 + h * dX[i]) / 6.0;
}

/// Return the derivatives (dx/dt dy/dt dz/dt)
void Lorenz::lorzrk(const double X[3], double deriv[3])
{
    deriv[0] = m_sigma * (X[1] - X[0]);
    deriv[1] = m_lambda * X[0] - X[1] - X[0] * X[2];
    deriv[2] = X[0] * X[1] - m_beta * X[2];
}

void Lorenz::output() const
{
    std::cout<<m_x[0]<<" "<<m_x[1]<<" "<<m_x[2]<<std::endl;
}
