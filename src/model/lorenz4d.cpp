#include "lorenz4d.h"

using namespace ldas;
Lorenz4d::Lorenz4d():m_F(10),m_h(0.01),m_nvar(40)
{
    //ctor
}

Lorenz4d::Lorenz4d(const double F, const double h, const unsigned int n)
    :m_F(F),m_h(h),m_nvar(n)
{

}

Lorenz4d::~Lorenz4d()
{
    //dtor
}

Lorenz4d::Lorenz4d(const Lorenz4d& other)
{
    //copy ctor
    m_h=0.01;
    m_F=10;
    m_X=other.parameter();
}

Lorenz4d& Lorenz4d::operator=(const Lorenz4d& other)
{
    if (this == &other) return *this; // handle self assignment
    //assignment operator
    m_h=0.01;
    m_F=10;
    m_X=other.parameter();
    return *this;
}

void Lorenz4d::parameter(const mat& X)
{
	m_X=X;
}
mat Lorenz4d::parameter() const
{
    return m_X;
}
mat Lorenz4d::result() const
{
    return m_X;
}

void Lorenz4d::run()
{
    rk4(m_h,m_F);
}

void Lorenz4d::rk4(const double h,const double F)
{
    mat d1=zeros<mat>(m_nvar,1);
    mat d2=zeros<mat>(m_nvar,1);
    mat d3=zeros<mat>(m_nvar,1);
    mat d4=zeros<mat>(m_nvar,1);

    lorzrk(m_X,F,d1);
    lorzrk(m_X+d1*0.5*h,F,d2);
    lorzrk(m_X+0.5*h*d2,F,d3);
    lorzrk(m_X+h*d3,F,d4);

    m_X+=1/6.0*h*(d1+2*d2+2*d3+d4);
}

/// Return the derivatives (dx/dt dy/dt dz/dt)
void Lorenz4d::lorzrk(const mat& X, const double F, mat& deriv)
{
	//int J=40;
	//deriv.zeros(m_nvar,1);
	deriv(0)=(X(1)-X(m_nvar-2))*X(m_nvar-1)-X(0);
	deriv(1)=(X(2)-X(m_nvar-1))*X(0)-X(1);
	deriv(m_nvar-1)=(X(0)-X(m_nvar-3))*X(m_nvar-2)-X(m_nvar-1);
	for(int i=2;i<m_nvar-1;i++)
	{
		deriv(i)=(X(i+1)-X(i-2))*X(i-1)-X(i);
	}
	deriv+=F;
}

void Lorenz4d::output() const
{
    std::cout<<trans(m_X);
}
