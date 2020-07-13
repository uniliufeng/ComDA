#include "SIR.h"

using namespace ldas;
SIR::SIR()
{
    //ctor
}

SIR::~SIR()
{
    //dtor
}

SIR::SIR(const SIR& other)
{
    //copy ctor
}

SIR& SIR::operator=(const SIR& other)
{
    if (this == &other) return *this; // handle self assignment

    return *this;
}

void SIR::parameter(const double beta, double gamma)
{
	m_beta = beta;
	m_gamma = gamma;
}

void SIR::init(const double N, double S, double I, double R)
{
	m_N = N;
	m_I = I;
	m_R = R;
	if (S==0)
	{m_S = m_N - m_I - m_R;}
	else
	{m_S = S;}
	
	
	m_state[0] = m_S;
	m_state[1] = m_I;
	m_state[2] = m_R;
}

void SIR::update(const double S, double I, double R)
{
	m_I = I;
	m_R = R;
	m_S = S;
	if (S==0)
	{m_S = m_N - m_I - m_R;}
	else
	{m_S = S;}
	
	m_state[0] = m_S;
	m_state[1] = m_I;
	m_state[2] = m_R;
}

void SIR::updateI(const double I)
{
	m_I = I;
	//m_S = m_N - m_I - m_E - m_R;
	
	//m_state[0] = m_S;
	m_state[1] = m_I;
	//m_state[2] = m_R;
}

void SIR::run()
{
	double rS = m_S;
	double rI = m_I;
	double rR = m_R;
	m_S = rS-m_beta*rI*rS/m_N;
	m_I = rI+m_beta*rI*rS/m_N;
	m_R = rR+m_gamma*rI;
	
	m_state[0] = m_S;
	m_state[1] = m_I;
	m_state[2] = m_R;
}

void SIR::output() const
{
    std::cout<<m_S<<" "<<m_I<<" "<<m_R<<std::endl;
}

const double* SIR::result() const
{
    return m_state;
}
